#include <fstream>
#include <iomanip>

#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>
#include <AMReX_PlotFileUtil.H>

#include <AMReX_BLFort.H>
#include <AMReX_Gpu.H>

#include <cvode/cvode.h>               /* prototypes for CVODE fcts., consts. */
#include <sunlinsol/sunlinsol_spgmr.h> /* access to SPGMR SUNLinearSolver     */
#include <cvode/cvode_spils.h>         /* access to CVSpils interface */
#include <cvode/cvode_diag.h>          /* access to CVDiag interface */
#include <sundials/sundials_types.h>   /* definition of type realtype */
#include <sundials/sundials_math.h>    /* definition of ABS and EXP   */

#include <nvector/nvector_cuda.h>
#include <nvector/nvector_serial.h>
#include "myfunc_F.H"

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

    // What time is it now?  We'll use this to compute total run time.
    Real strt_time = ParallelDescriptor::second();

    std::cout << std::setprecision(15);

    int n_cell, max_grid_size;
    int cvode_meth, cvode_itmeth, write_plotfile;
    bool do_tiling;

  realtype reltol, abstol, t, tout, umax;
  N_Vector u;
  //  SUNLinearSolver LS;
  void *cvode_mem;
  int iout, flag;
  long int nst;

  u = NULL;
  //  LS = NULL;
  cvode_mem = NULL;

  reltol = 1e-6;  /* Set the tolerances */
  abstol = 1e-10;

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

      // Select CVODE solve method.
      //   1 = Adams (for non-stiff problems)
      //   2 = BDF (for stiff problems)
      pp.get("cvode_meth",cvode_meth);
      // Select CVODE solver iteration method.
      //   1 = Functional iteration
      //   2 = Newton iteration
      pp.get("cvode_itmeth",cvode_itmeth);

      pp.get("write_plotfile",write_plotfile);
      pp.get("do_tiling",do_tiling);
    }

    if (cvode_meth < 1)
      amrex::Abort("Unknown cvode_meth");
    if (cvode_itmeth < 1)
      amrex::Abort("Unknown cvode_itmeth");

    amrex::Print() << "This is AMReX version " << amrex::Version() << std::endl;
    amrex::Print() << "Problem domain size: nx = ny = nz = " << n_cell << std::endl;
    amrex::Print() << "Max grid size: " << max_grid_size << std::endl;
    amrex::Print() << "CVODE method: ";
    if (cvode_meth == 1) {
      amrex::Print() << "Adams (non-stiff)";
    } else if (cvode_meth == 2) {
        amrex::Print() << "BDF (stiff)";
    }
    amrex::Print() << std::endl;
    amrex::Print() << "CVODE iteration method: ";
    if (cvode_itmeth == 1) {
      amrex::Print() << "Functional";
    } else if (cvode_itmeth == 2) {
        amrex::Print() << "Newton";
    }
    amrex::Print() << std::endl;

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

    // Create MultiFab with no ghost cells.
    MultiFab mf(ba, dm, Ncomp, 0);

    
    mf.setVal(0.0);

    //    MultiFab vode_aux_vars(ba, dm, 4*Ncomp, 0);
      /*      vode_aux_vars.setVal(3.255559960937500E+04,0);
      vode_aux_vars.setVal(1.076699972152710E+00,1);
      vode_aux_vars.setVal(2.119999946752000E+12,2);
      vode_aux_vars.setVal(1/(1+1.635780036449432E-01),3);*/
    //////////////////////////////////////////////////////////////////////
    // Allocate data
    //////////////////////////////////////////////////////////////////////
    
    fort_init_allocations();
    fort_init_tables_eos_params();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for ( MFIter mfi(mf, do_tiling); mfi.isValid(); ++mfi )
    {
      t=0;
      tout=2;
      tout=8.839029760565609E-06;

      Real* dptr;

      const Box& tbx = mfi.tilebox();
      amrex::IntVect tile_size = tbx.size();
      const int* hi = tbx.hiVect();
      const int* lo = tbx.loVect();
      int long neq=(tile_size[0])*(tile_size[1])*(tile_size[2]);

      if(neq>1)
	amrex::Print()<<"Integrating a box with "<<neq<<" cels"<<std::endl;

      /* Create a CUDA vector with initial values */
      u = N_VNew_Cuda(neq);  /* Allocate u vector */
      if(check_retval((void*)u, "N_VNew_Cuda", 0)) return(1);

      FSetInternalEnergy_mfab(mf[mfi].dataPtr(),
        tbx.loVect(),
	    tbx.hiVect());  /* Initialize u vector */

      dptr=N_VGetHostArrayPointer_Cuda(u);
      amrex::Cuda::Device::synchronize();
      AMREX_GPU_ERROR_CHECK();
      mf[mfi].copyToMem(tbx,0,1,dptr);
      N_VCopyToDevice_Cuda(u);


      /* Call CVodeCreate to create the solver memory and specify the 
       * Backward Differentiation Formula and the use of a Newton iteration */
      #ifdef CV_NEWTON
      cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
      #else
      cvode_mem = CVodeCreate(CV_BDF);
      #endif
      if(check_retval((void *)cvode_mem, "CVodeCreate", 0)) return(1);

      amrex::Cuda::Device::synchronize();
      AMREX_GPU_ERROR_CHECK();
      /* Call CVodeInit to initialize the integrator memory and specify the
       * user's right hand side function in u'=f(t,u), the initial time T0, and
       * the initial dependent variable vector u. */
      flag = CVodeInit(cvode_mem, f, t, u);
      if(check_retval(&flag, "CVodeInit", 1)) return(1);

      /* Call CVodeSStolerances to specify the scalar relative tolerance
       * and scalar absolute tolerance */
      N_Vector abstol_vec = N_VNew_Cuda(neq);
      N_VConst(abstol,abstol_vec);
      flag = CVodeSVtolerances(cvode_mem, reltol, abstol_vec);
      //flag = CVodeSStolerances(cvode_mem, reltol, abstol);
      if (check_retval(&flag, "CVodeSVtolerances", 1)) return(1);

      /* Create SPGMR solver structure without preconditioning
       * and the maximum Krylov dimension maxl */
      //      LS = SUNSPGMR(u, PREC_NONE, 0);
      //      if(check_flag(&flag, "SUNSPGMR", 1)) return(1);

      /* Set CVSpils linear solver to LS */
      //      flag = CVSpilsSetLinearSolver(cvode_mem, LS);
      //      if(check_flag(&flag, "CVSpilsSetLinearSolver", 1)) return(1);

      flag = CVDiag(cvode_mem);

      /*Use N_Vector to create userdata, in order to allocate data on device*/
      N_Vector Data = N_VNew_Cuda(4*neq);  // Allocate u vector 
      N_VConst(0.0,Data);
      amrex::Cuda::Device::synchronize();
      AMREX_GPU_ERROR_CHECK();
      double* rparh=N_VGetHostArrayPointer_Cuda(Data);
      for(int i=0;i<neq;i++)
	{
	  rparh[4*i+0]= 3.255559960937500E+04;   //rpar(1)=T_vode
	  rparh[4*i+1]= 1.076699972152710E+00;//    rpar(2)=ne_vode
	  rparh[4*i+2]=  2.119999946752000E+12; //    rpar(3)=rho_vode
	  rparh[4*i+3]=1/(1.635780036449432E-01)-1;    //    rpar(4)=z_vode

	}
      N_VCopyToDevice_Cuda(Data);
      //      amrex::Device::synchronize();
      //      CudaErrorCheck();
      /////      CVodeSetUserData(cvode_mem, N_VGetHostArrayPointer_Cuda(Data));
      CVodeSetUserData(cvode_mem, &Data);
      //      CVodeSetUserData(cvode_mem, N_VGetDeviceArrayPointer_Cuda(Data));
      /*      double* dptr_data=new double[4];
      for(int i=0;i<4;i++)
      dptr_data[i]=0.0;
      CVodeSetUserData(cvode_mem, dptr_data);*/

      /* Call CVode */
      /*      flag = CVode(cvode_mem, tout, u, &t, CV_NORMAL);
	      if(check_flag(&flag, "CVode", 1)) break;*/
      for(iout=1, tout=8.839029760565609E-06  ; iout <= 1; iout++, tout += 8.839029760565609E-06) {
      flag = CVode(cvode_mem, tout, u, &t, CV_NORMAL);
      umax = N_VMaxNorm(u);
      N_VCopyFromDevice_Cuda(Data);
      /*
      double min_ne=10;
      double max_ne=-1;
      for(int i=0;i<neq;i++)
	{
	min_ne=min(min_ne,rparh[4*i+1]);
	max_ne=max(max_ne,rparh[4*i+1]);
      }
      amrex::Print() << "min is" << min_ne<< "\tmax is "<< max_ne << std::endl;*/
      flag = CVodeGetNumSteps(cvode_mem, &nst);
      check_retval(&flag, "CVodeGetNumSteps", 1);
      PrintOutput(tout, umax, nst);
      }

      amrex::Cuda::Device::synchronize();
      AMREX_GPU_ERROR_CHECK();

      N_VCopyFromDevice_Cuda(u);
      /*      N_VCopyFromDevice_Cuda(Data);
	      fprintf(stdout,"\nFinal rparh[0]=%g \n\n",rparh[0]);*/
      /////
      mf[mfi].copyFromMem(tbx,0,1,dptr);
      /////      CVodeSetUserData(cvode_mem, NULL);
      PrintFinalStats(cvode_mem);

      amrex::Cuda::Device::synchronize();
      AMREX_GPU_ERROR_CHECK();

      N_VDestroy(u);          /* Free the u vector */
      /*      delete(dptr_data);
      dptr_data=NULL;*/
      N_VDestroy(Data);          /* Free the userdata vector */
      CVodeFree(&cvode_mem);  /* Free the integrator memory */
    
    }

    //////////////////////////////////////////////////////////////////////
    // Deallocate data
    //////////////////////////////////////////////////////////////////////
    
    fort_fin_allocations();

    amrex::Print()<<"Maximum of repacked final solution: "<<mf.max(0,0,0)<<std::endl;
    
    if (write_plotfile)
    {
      amrex::WriteSingleLevelPlotfile("PLT_OUTPUT",
                                      mf,
                                      {"y1"},
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
    
    amrex::Finalize();
    return 0;
}

__global__ void f_rhs_test(Real t,double* u_ptr,Real* udot_ptr, Real* rpar, int neq)
{
  /*
1.635780036449432E-01 a
8.839029760565609E-06 dt
2.119999946752000E+12 rho
3.255559960937500E+04 T
1.076699972152710E+00 ne
6.226414794921875E+02 e */
  /*  double rpar2[4];
    rpar2[0]= 3.255559960937500E+04;   //rpar(1)=T_vode
  rpar2[1]= 1.076699972152710E+00;//    rpar(2)=ne_vode
  rpar2[2]=  2.119999946752000E+12; //    rpar(3)=rho_vode
  rpar2[3]=1/(1+1.635780036449432E-01);    //    rpar(4)=z_vode*/
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
  //  if(tid==0)
    if(tid<neq)
      RhsFn(t,u_ptr+tid,udot_ptr+tid,rpar+4*tid,1);
    //    udot_ptr[tid]=(neq-tid)*t
  
    //*********************************
    /*    if(tid<neq)
	  RhsFn(t,u_ptr+tid,udot_ptr+tid,rpar+4*tid,1);*/
    //rpar[4*tid+1]=tid;
    /**********************************
    if(tid<neq)
      udot_ptr[tid]=2.0*t;*/

  /* Either way to setup IC seems to work
    RhsFn(t,u_ptr+i,udot_ptr+i,rpar2,1);*/
}

static int f(realtype t, N_Vector u, N_Vector udot, void *user_data)
{
  Real* udot_ptr=N_VGetDeviceArrayPointer_Cuda(udot);
  Real* u_ptr=N_VGetDeviceArrayPointer_Cuda(u);
  int neq=N_VGetLength_Cuda(udot);
  double*  rpar=N_VGetDeviceArrayPointer_Cuda(*(static_cast<N_Vector*>(user_data)));
  
  /*  N_VCopyFromDevice_Cuda(*(static_cast<N_Vector*>(user_data)));  
  double*  rparh=N_VGetHostArrayPointer_Cuda(*(static_cast<N_Vector*>(user_data)));
   
   fprintf(stdout,"\nt=%g \n\n",t);
   fprintf(stdout,"\nrparh[0]=%g \n\n",rparh[0]);
  fprintf(stdout,"\nrparh[1]=%g \n\n",rparh[1]);
  fprintf(stdout,"\nrparh[2]=%g \n\n",rparh[2]);
  fprintf(stdout,"\nrparh[3]=%g \n\n",rparh[3]);*/
  int numThreads = std::min(32, neq);
  int numBlocks = static_cast<int>(ceil(((double) neq)/((double) numThreads)));
  /*
 ////////////////////////////  fprintf(stdout,"\n castro <<<%d,%d>>> \n\n",numBlocks, numThreads);

    unsigned block = 256;
    unsigned grid = (int) ceil((float)neq / block);
 ////////////////////////////fprintf(stdout,"\n cvode <<<%d,%d>>> \n\n",grid, block);

 int blockSize, gridSize;
 
 // Number of threads in each thread block
 blockSize = 1024;
 
 // Number of thread blocks in grid
 gridSize = (int)ceil((float)neq/blockSize);
 ////////////////////////////fprintf(stdout,"\n olcf <<<%d,%d>>> \n\n",gridSize, blockSize);
 */
 f_rhs_test<<<numBlocks,numThreads>>>(t,u_ptr,udot_ptr, rpar, neq);

 amrex::Cuda::Device::synchronize();
 AMREX_GPU_ERROR_CHECK();

/*      N_VCopyFromDevice_Cuda(*(static_cast<N_Vector*>(user_data)));
      fprintf(stdout,"\nafter rparh[0]=%g \n\n",rparh[0]);
  fprintf(stdout,"\nafter rparh[1]=%g \n\n",rparh[1]);
  fprintf(stdout,"\nafter rparh[2]=%g \n\n",rparh[2]);
  fprintf(stdout,"\nafter rparh[3]=%g \n\n",rparh[3]);
  fprintf(stdout,"\nafter last rparh[4*(neq-1)+1]=%g \n\n",rparh[4*(neq-1)+1]);*/
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
