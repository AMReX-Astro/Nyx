#include <fstream>
#include <iomanip>

#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>
#include <AMReX_PlotFileUtil.H>
//#if !defined(BL_NO_FORT)
//#include <AMReX_BaseFab_f.H>
//#endif

#include <AMReX_BLFort.H>

#include <cvode/cvode.h>               /* prototypes for CVODE fcts., consts. */
#include <cvode/cvode_diag.h>          /* access to CVDiag interface */
#include <sundials/sundials_types.h>   /* definition of type realtype */
#include <sundials/sundials_math.h>    /* definition of ABS and EXP   */

#include <nvector/nvector_serial.h>
#ifdef _OPENMP
#include <nvector/nvector_openmp.h>
#endif
#ifdef AMREX_USE_CUDA
#include <nvector/nvector_cuda.h>
#endif
#include "myfunc_F.H"

  int integrate_state_vec(amrex::MultiFab &state,   amrex::MultiFab &diag_eos, const amrex::Real& a, const amrex::Real& delta_time);
/* Functions Called by the Solver */
static int f(realtype t, N_Vector u, N_Vector udot, void *user_data);

static void PrintFinalStats(void *cvode_mem);

static void PrintOutput(realtype t, realtype umax, long int nst);

/* Private function to check function return values */
static int check_retval(void *flagvalue, const char *funcname, int opt);

using namespace amrex;

/*
#pragma gpu
void device_ptr_wrapper(N_Vector u, Real* dptr)
      dptr=N_VGetDeviceArrayPointer_Cuda(u);
*/

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    {

    // What time is it now?  We'll use this to compute total run time.
    Real strt_time = ParallelDescriptor::second();

    std::cout << std::setprecision(15);

    bool test_ic;
    int n_cell, max_grid_size;
    int write_plotfile;
    bool do_tiling;

    // inputs parameters
    {
      // ParmParse is way of reading inputs from the inputs file
      ParmParse pp;

      // We need to get n_cell from the inputs file - this is the number of
      // cells on each side of a square (or cubic) domain.
      pp.get("n_cell",n_cell);

      // Default nsteps to 0, allow us to set it to something else in the
      // inputs file
      pp.get("max_grid_size",max_grid_size);

      pp.get("write_plotfile",write_plotfile);
      pp.get("do_tiling",do_tiling);
    }

    amrex::Print() << "This is AMReX version " << amrex::Version() << std::endl;
    amrex::Print() << "Problem domain size: nx = ny = nz = " << n_cell << std::endl;
    amrex::Print() << "Max grid size: " << max_grid_size << std::endl;

    // make BoxArray and Geometry
    BoxArray ba;
    Geometry geom;
    {
      IntVect dom_lo(IntVect(D_DECL(0,0,0)));
      IntVect dom_hi(IntVect(D_DECL(n_cell-1, n_cell-1, n_cell-1)));
      Box domain(dom_lo, dom_hi);

      // Initialize the boxarray "ba" from the single box "bx"
      ba.define(domain);

      // Break up boxarray "ba" into chunks no larger than "max_grid_size"
      // along a direction
      ba.maxSize(max_grid_size);

      // This defines the physical size of the box.  Right now the box is
      // [-1,1] in each direction.
      RealBox real_box;
      for (int n = 0; n < BL_SPACEDIM; n++) {
        real_box.setLo(n,-1.0);
        real_box.setHi(n, 1.0);
      }

      // This sets the boundary conditions to be doubly or triply periodic
      int is_periodic[BL_SPACEDIM];
      for (int i = 0; i < BL_SPACEDIM; i++) {
        is_periodic[i] = 1;
      }

      // This defines a Geometry object
      geom.define(domain,&real_box,CoordSys::cartesian,is_periodic);
    }

    // Ncomp = number of components for each array
    int Ncomp  = 1;

    // time = starting time in the simulation
    Real time = 0.0;

    DistributionMapping dm(ba);

    //If ghost cells are needed, growntilebox must be used in integrator
    // Create MultiFab with no ghost cells.
    MultiFab S_old(ba, dm, 7, 0);

    // Create MultiFab with no ghost cells.
    MultiFab D_old(ba, dm, 4, 0);

    Real a = 1.635780036449432E-01;
    S_old.setVal(2.119999946752000E+12,0,1,0); //    rpar(3)=rho_vode
    S_old.setVal(2.119999946752000E+12*6.226414794921875E+02,5,1,0); //    rho*e
    S_old.setVal(2.119999946752000E+12*(6.226414794921875E+02+2E+04),4,1,0); //    rho*E
    D_old.setVal(3.255559960937500E+04 ,0,1,0); //    rpar(1)=T_vode
    D_old.setVal(1.076699972152710E+00 ,1,1,0); //    rpar(2)=ne_vode
    Real half_dt=1e-5;
    //////////////////////////////////////////////////////////////////////
    // Allocate data
    //////////////////////////////////////////////////////////////////////
    
    AMREX_GPU_ERROR_CHECK();
    fort_alloc_cuda_managed();
    fort_init_tables_eos_params();
    AMREX_GPU_ERROR_CHECK();
    int ierr=integrate_state_vec(S_old,       D_old,       a, half_dt);
    //////////////////////////////////////////////////////////////////////
    // Deallocate data
    //////////////////////////////////////////////////////////////////////
    
    fort_dealloc_cuda_managed();

    amrex::Print()<<"Maximum of repacked final solution: "<<S_old.max(5,0,0)<<std::endl;
    
    if (write_plotfile)
    {
      amrex::Vector<std::string, std::allocator<std::string>> varnames({"rho","xmom","ymom","zmom","rho_E","rho_e","6"});
      amrex::Vector<std::string, std::allocator<std::string>> varnames_diag({"Temp","Ne","3","4"});
      amrex::WriteSingleLevelPlotfile("PLT_OUTPUT",
                                      S_old,
				      varnames,
                                      geom,
                                      time,
                                      0);
      amrex::WriteSingleLevelPlotfile("PLT_DIAG_OUTPUT",
                                      D_old,
				      varnames_diag,
                                      geom,
                                      time,
                                      0);
    }

    // Call the timer again and compute the maximum difference between the start time and stop time
    //   over all processors
    Real stop_time = ParallelDescriptor::second() - strt_time;
    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    ParallelDescriptor::ReduceRealMax(stop_time,IOProc);

    // Tell the I/O Processor to write out the "run time"
    amrex::Print() << "Run time = " << stop_time << std::endl;
    }
    amrex::Finalize();
    return 0;
}

int integrate_state_vec
  (amrex::MultiFab &S_old,
   amrex::MultiFab &D_old,
   const Real& a, const Real& delta_time)
{
    // time = starting time in the simulation
  realtype reltol, abstol, tout, umax;
  N_Vector u;
  //  SUNLinearSolver LS;
  int iout, flag;
  bool do_tiling=false;    

  reltol = 1e-4;  /* Set the tolerances */
  abstol = 1e-4;

  /*
  int Temp_loc=Temp_comp;
  int Ne_loc=Ne_comp;
  int Eint_loc=Eint;
  int Eden_loc=Eden;
  int Density_loc=Density;
  int one_in = 1;
  */
  int Temp_loc=0;
  int Ne_loc=1;
  int Eint_loc=5;
  int Eden_loc=4;
  int Density_loc=0;
  int one_in = 1;

  fort_ode_eos_setup(a,delta_time);

// This pragma version would split up the mfiter loops
// and assign an OMP thread per tile, in which case the serial
// NVector should be used
  //#ifdef _OPENMP
  //#pragma omp parallel if (Gpu::notInLaunchRegion())
  //#endif
#ifdef _OPENMP
  for ( MFIter mfi(S_old, false); mfi.isValid(); ++mfi )
#else
  for ( MFIter mfi(S_old, TilingIfNotGPU()); mfi.isValid(); ++mfi )
#endif
    {

      double* dptr;
      //check that copy contructor vs create constructor works??
      const Box& tbx = mfi.tilebox();
      //Array4<Real> const& state = S_old.array(mfi);
      //Array4<Real> const& diag_eos = D_old.array(mfi);
      const auto len = amrex::length(tbx);  // length of box
      const auto lo  = amrex::lbound(tbx);  // lower bound of box
      const auto state = (S_old[mfi]).view(lo);  // a view starting from lo
      const auto diag_eos = (D_old[mfi]).view(lo);  // a view starting from lo

      Array4<Real> const& state4 = S_old.array(mfi);
      Array4<Real> const& diag_eos4 = D_old.array(mfi);

      int width = len.x; //trying to run all i
      width = len.x*len.y*len.z;

      /*      for       (int k = 0; k < len.z; ++k) {
	  // We know this is safe for simd on cpu.  So let's give compiler some help.
	//	AMREX_PRAGMA_SIMD
	for     (int j = 0; j < len.y; ++j) {
	for (int i = 0; i < len.x; i+=width) {*/
				N_Vector u;
				//  SUNLinearSolver LS;
				void *cvode_mem;
				realtype t=0.0;
				
				u = NULL;
				//  LS = NULL;
				cvode_mem = NULL;
				//				long int neq = len.x;
				long int neq = width;
				int loop = 1;

#ifdef AMREX_USE_CUDA
				cudaStream_t currentStream = amrex::Cuda::Device::cudaStream();
				u = N_VNew_Cuda(neq);  /* Allocate u vector */
				N_Vector e_orig = N_VNew_Cuda(neq);  /* Allocate u vector */
				N_VSetCudaStream_Cuda(e_orig, &currentStream);
				double* eptr=N_VGetDeviceArrayPointer_Cuda(e_orig);
				N_VSetCudaStream_Cuda(u, &currentStream);
				dptr=N_VGetDeviceArrayPointer_Cuda(u);

				N_Vector Data = N_VNew_Cuda(4*neq);  // Allocate u vector 
				N_VSetCudaStream_Cuda(Data, &currentStream);
				N_VConst(0.0,Data);
				double* rparh=N_VGetDeviceArrayPointer_Cuda(Data);
				N_Vector abstol_vec = N_VNew_Cuda(neq);
				N_VSetCudaStream_Cuda(abstol_vec,&currentStream);
#else
#ifdef _OPENMP
				int nthreads=omp_get_max_threads();
				u = N_VNew_OpenMP(neq,nthreads);  /* Allocate u vector */
				N_Vector e_orig = N_VNew_OpenMP(neq,nthreads);  /* Allocate u vector */
				double* eptr=N_VGetArrayPointer_Serial(e_orig);
				dptr=N_VGetArrayPointer_OpenMP(u);

				N_Vector Data = N_VNew_OpenMP(4*neq,nthreads);  // Allocate u vector 
				N_VConst(0.0,Data);
				double* rparh=N_VGetArrayPointer_OpenMP(Data);
				N_Vector abstol_vec = N_VNew_OpenMP(neq,nthreads);
#else
				u = N_VNew_Serial(neq);  /* Allocate u vector */
				N_Vector e_orig = N_VNew_Serial(neq);  /* Allocate u vector */
				double* eptr=N_VGetArrayPointer_Serial(e_orig);
				dptr=N_VGetArrayPointer_Serial(u);

				N_Vector Data = N_VNew_Serial(4*neq);  // Allocate u vector 
				N_VConst(0.0,Data);
				double* rparh=N_VGetArrayPointer_Serial(Data);
				N_Vector abstol_vec = N_VNew_Serial(neq);
#endif
#endif	

#ifdef _OPENMP
      const Dim3 hi = amrex::ubound(tbx);

#pragma omp parallel for collapse(3)
      for (int k = lo.z; k <= hi.z; ++k) {
	for (int j = lo.y; j <= hi.y; ++j) {
	    for (int i = lo.x; i <= hi.x; ++i) {
	      fort_ode_eos_setup(a,delta_time);
#else
				AMREX_PARALLEL_FOR_3D ( tbx, i,j,k,
				{				  
#endif
				  // This index is a mix of the array4 and the lo view
				  int idx = i+j*len.x+k*len.x*len.y-(lo.x+lo.y*len.x+lo.z*len.x*len.y);
				  dptr[idx*loop]=state4(i,j,k,Eint_loc)/state4(i,j,k,Density_loc);
				  eptr[idx*loop]=state4(i,j,k,Eint_loc)/state4(i,j,k,Density_loc);
				  //Consider pre-subtracting original internal energy contribution
				  if(false)
				    {
				      state4(i,j,k,Eint_loc)  -= state4(i,j,k,Density_loc) * (dptr[idx*loop]);
				      state4(i,j,k,Eden_loc)  -= state4(i,j,k,Density_loc) * (dptr[idx*loop]);
				    }
				  rparh[4*idx*loop+0]= diag_eos4(i,j,k,Temp_loc);   //rpar(1)=T_vode
				  rparh[4*idx*loop+1]= diag_eos4(i,j,k,Ne_loc);//    rpar(2)=ne_vode
				  rparh[4*idx*loop+2]= state4(i,j,k,Density_loc); //    rpar(3)=rho_vode
				  rparh[4*idx*loop+3]=1/a-1;    //    rpar(4)=z_vode
				  //				}
#ifdef _OPENMP
				}
				}
				}
#pragma omp barrier
#else
				});
#endif

#ifdef CV_NEWTON
				cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
#else
				cvode_mem = CVodeCreate(CV_BDF);
#endif
				flag = CVodeInit(cvode_mem, f, t, u);

				N_VScale(abstol,u,abstol_vec);
				flag = CVodeSVtolerances(cvode_mem, reltol, abstol_vec);

				//				flag = CVodeSStolerances(cvode_mem, reltol, dptr[0]*abstol);
				flag = CVDiag(cvode_mem);

				CVodeSetMaxNumSteps(cvode_mem,2000);

				N_Vector constrain=N_VClone(u);
				N_VConst(2,constrain);	      
				flag =CVodeSetConstraints(cvode_mem,constrain);
				
				CVodeSetUserData(cvode_mem, &Data);
				//				CVodeSetMaxStep(cvode_mem, delta_time/10);
				BL_PROFILE_VAR("Nyx::strang_second_cvode",cvode_timer2);
				flag = CVode(cvode_mem, delta_time, u, &t, CV_NORMAL);
				amrex::Cuda::Device::streamSynchronize();
				BL_PROFILE_VAR_STOP(cvode_timer2);

#ifndef NDEBUG
				PrintFinalStats(cvode_mem);
#endif

#ifdef _OPENMP
#pragma omp parallel for collapse(3)
      for (int k = lo.z; k <= hi.z; ++k) {
	for (int j = lo.y; j <= hi.y; ++j) {
	    for (int i = lo.x; i <= hi.x; ++i) {
#else
				AMREX_PARALLEL_FOR_3D ( tbx, i,j,k,
				{				  
#endif
				  int  idx= i+j*len.x+k*len.x*len.y-(lo.x+lo.y*len.x+lo.z*len.x*len.y);
				  fort_ode_eos_finalize(&(dptr[idx*loop]), &(rparh[4*idx*loop]), one_in);
				  diag_eos4(i,j,k,Temp_loc)=rparh[4*idx*loop+0];   //rpar(1)=T_vode
				  diag_eos4(i,j,k,Ne_loc)=rparh[4*idx*loop+1];//    rpar(2)=ne_vode
				
				  state4(i,j,k,Eint_loc)  += state4(i,j,k,Density_loc) * (dptr[idx*loop]-eptr[idx]);
				  state4(i,j,k,Eden_loc)  += state4(i,j,k,Density_loc) * (dptr[idx*loop]-eptr[idx]);
#ifdef _OPENMP
				}
				}
				}
#pragma omp barrier
#else
				});
#endif

				N_VDestroy(u);          /* Free the u vector */
				N_VDestroy(abstol_vec);          /* Free the u vector */
				N_VDestroy(Data);          /* Free the userdata vector */
				CVodeFree(&cvode_mem);  /* Free the integrator memory */

    }

  #ifdef NDEBUG
  #ifndef NDEBUG
        if (S_old.contains_nan())
            amrex::Abort("state has NaNs after the second strang call");
  #endif
#endif

    return 0;
}

static int f(realtype t, N_Vector u, N_Vector udot, void *user_data)
{
#ifdef AMREX_USE_CUDA
  Real* udot_ptr=N_VGetDeviceArrayPointer_Cuda(udot);
  Real* u_ptr=N_VGetDeviceArrayPointer_Cuda(u);
  int neq=N_VGetLength_Cuda(udot);
  double*  rpar=N_VGetDeviceArrayPointer_Cuda(*(static_cast<N_Vector*>(user_data)));
#else
  Real* udot_ptr=N_VGetArrayPointer_Serial(udot);
  Real* u_ptr=N_VGetArrayPointer_Serial(u);
  int neq=N_VGetLength_Serial(udot);
  double*  rpar=N_VGetArrayPointer_Serial(*(static_cast<N_Vector*>(user_data)));
#endif
  
  int numThreads = std::min(32, neq);
  int numBlocks = static_cast<int>(ceil(((double) neq)/((double) numThreads)));
  //  cudaStream_t currentStream = amrex::Cuda::Device::cudaStream();
#ifdef _OPENMP
  #pragma omp parallel for
  for(int idx=0;idx<neq;idx++)
      RhsFnReal(t,u_ptr+idx,udot_ptr+idx, rpar+4*idx, 1);
#else
  AMREX_LAUNCH_DEVICE_LAMBDA ( neq, idx, {
      RhsFnReal(t,u_ptr+idx,udot_ptr+idx, rpar+4*idx, 1);
  });
#endif
  amrex::Cuda::Device::streamSynchronize();//cudaStreamSynchronize(currentStream);
  AMREX_GPU_ERROR_CHECK();

  return 0;
}

static void PrintOutput(realtype t, realtype umax, long int nst)
{
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("At t = %4.2Lf   max.norm(u) =%14.16Le   nst = %4ld\n", t, umax, nst);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("At t = %4.2f   max.norm(u) =%14.16e   nst = %4ld\n", t, umax, nst);
#else
  printf("At t = %4.2f   max.norm(u) =%14.16e   nst = %4ld\n", t, umax, nst);
#endif

  return;
}

/* Get and print some final statistics */

static void PrintFinalStats(void *cvode_mem)
{
  long lenrw, leniw ;
  long lenrwLS, leniwLS;
  long int nst, nfe, nsetups, nni, ncfn, netf;
  long int nli, npe, nps, ncfl, nfeLS;
  int retval;

  retval = CVodeGetWorkSpace(cvode_mem, &lenrw, &leniw);
  check_retval(&retval, "CVodeGetWorkSpace", 1);
  retval = CVodeGetNumSteps(cvode_mem, &nst);
  check_retval(&retval, "CVodeGetNumSteps", 1);
  retval = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  check_retval(&retval, "CVodeGetNumRhsEvals", 1);
  retval = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  check_retval(&retval, "CVodeGetNumLinSolvSetups", 1);
  retval = CVodeGetNumErrTestFails(cvode_mem, &netf);
  check_retval(&retval, "CVodeGetNumErrTestFails", 1);
  retval = CVDiagGetNumRhsEvals(cvode_mem, &nfeLS);

  printf("\nFinal Statistics.. \n\n");
  printf("lenrw   = %5ld     leniw   = %5ld\n"  , lenrw, leniw);
  printf("nst     = %5ld\n"                     , nst);
  printf("nfe     = %5ld     nfeLS   = %5ld\n"  , nfe, nfeLS);
  printf("nsetups = %5ld     netf    = %5ld\n"  , nsetups, netf);

  return;
}

static int check_retval(void *flagvalue, const char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */

  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return(1); }

  /* Check if flag < 0 */

  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
              funcname, *errflag);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */

  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return(1); }

  return(0);
}
