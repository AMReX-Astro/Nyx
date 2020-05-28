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
#include <arkode/arkode_arkstep.h>      // prototypes for ARKStep fcts., consts.
#include <cvode/cvode_diag.h>          /* access to CVDiag interface */
#include <sunmatrix/sunmatrix_dense.h>     /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h>     /* access to dense SUNLinearSolver      */
#include <sunlinsol/sunlinsol_pcg.h>  /* access to PCG SUNLinearSolver        */
#include <sunlinsol/sunlinsol_spgmr.h>  /* access to SPGMR SUNLinearSolver             */
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

#include <iostream>

#ifdef AMREX_USE_CUDA
#include <cuda_profiler_api.h>
#endif

#include "SdcIntegrator.H"
#include "SparseGaussJordan.H"
#include "vode_system.H"
#include "RealVector.H"
#include "WallTimer.H"

enum State_Type_Index {
  Density = 0,
  Xmom    = 1,
  Ymom    = 2,
  Zmom    = 3,
  Eden    = 4,
  Eint    = 5

};

enum DiagEOS_Type_Index {
  Temp_comp = 0,
  Ne_comp   = 1

};

    int use_cvode=1;
    int use_cvdiag=1;
    bool do_tiling=false;

  int integrate_state_vec(amrex::MultiFab &state,   amrex::MultiFab &diag_eos, const amrex::Real& a, const amrex::Real& delta_time);
/* Functions Called by the Solver */
static int f(realtype t, N_Vector u, N_Vector udot, void *user_data);

static void PrintFinalStats(void *cvode_mem);

static void PrintOutput(realtype t, realtype umax, long int nst);

/* Private function to check function return values */
static int check_retval(void *flagvalue, const char *funcname, int opt);

using namespace amrex;

void* sunalloc(size_t mem_size)
{
  return (void*) The_Managed_Arena()->alloc(mem_size);
}

void sunfree(void* ptr)
{
  The_Managed_Arena()->free(ptr);
}

int integrate_state_vec_mfin(amrex::Array4<amrex::Real>const& state4,   amrex::Array4<amrex::Real>const& diag_eos4,const  amrex::Box& tbx,  const amrex::Real& a, const amrex::Real& delta_time);

int integrate_state_vec
  (amrex::MultiFab &S_old,
   amrex::MultiFab &D_old,
   const Real& a, const Real& delta_time)
{
    // time = starting time in the simulation

  amrex::Gpu::LaunchSafeGuard lsg(true);
  fort_ode_eos_setup(a,delta_time);

  //#ifdef _OPENMP
  //#pragma omp parallel if (Gpu::notInLaunchRegion())
  //#endif
#ifdef _OPENMP
  for ( MFIter mfi(S_old, do_tiling); mfi.isValid(); ++mfi )
#else
    //    for ( MFIter mfi(S_old, IntVect(AMREX_D_DECL(4,4,4))); mfi.isValid(); ++mfi )
    for ( MFIter mfi(S_old, do_tiling && TilingIfNotGPU()); mfi.isValid(); ++mfi )
#endif
    {

      //check that copy contructor vs create constructor works??
      const Box& tbx = mfi.tilebox();

      S_old[mfi].prefetchToDevice();
      D_old[mfi].prefetchToDevice();
      Array4<Real> const& state4 = S_old.array(mfi);
      Array4<Real> const& diag_eos4 = D_old.array(mfi);

      integrate_state_vec_mfin(state4,diag_eos4,tbx,a,delta_time);
}
      return 0;
}

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
    int test_match=1;
    bool increase_tol=0;

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
      pp.query("test_match",test_match);
      pp.query("use_cvode",use_cvode);
      pp.query("use_cvdiag",use_cvdiag);
      pp.query("increase_tol",increase_tol);
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

    Real a = 0.0688707121;

    if(test_match==0)
      {
	S_old.setVal(2.587879236e+10,0,1,0); //    rpar(3)=rho_vode
	S_old.setVal(2.920381766e+12,5,1,0); //    rho*e
      }
    else if(test_match==1)
      {
	S_old.setVal(1.83E12,5,1,0); //    rpar(3)=rho_vode
	S_old.setVal(1.98E10,0,1,0); //    rho*e
      }

    S_old.setVal(6.896610338e+12,4,1,0); //    rho*E
    S_old.setVal(0,6,1,0); //    rho*E
    D_old.setVal(11026.08482 ,0,1,0); //    rpar(1)=T_vode
    D_old.setVal(0.008728890037 ,1,1,0); //    rpar(2)=ne_vode
    Real half_dt=4.907133436e-06;
    //Real half_dt=1.0;
    //Real half_dt=1.6842553218299e-10;

	for (MFIter mfi(S_old,do_tiling); mfi.isValid(); ++mfi)
	  {
	    // Note that this "bx" includes the grow cells 
	    const Box& bx = mfi.tilebox();
	    ////	    IntVect test_point(AMREX_D_DECL(2,3,1));
	    IntVect test_point(AMREX_D_DECL(0,0,0));
	    if(bx.contains(test_point))
	      {
		const auto fab_state = S_old.array(mfi);
		if(test_match)
		  {
		    fab_state(test_point,5)=2.920381766e+12;//,5,1,0); //    rho*e
		    fab_state(test_point,6)=increase_tol;//,5,1,0); //    rho*e
		    fab_state(test_point,0)=2.587879236e+10;//,0,1,0); //    rpar(3)=rho_vode
		  }
		else
		  {
		    fab_state(test_point,5)=1.83E12;//,5,1,0); //    rho*e
		    fab_state(test_point,6)=increase_tol;//,5,1,0); //    rho*e
		    fab_state(test_point,0)=1.98E10;//,0,1,0); //    rpar(3)=rho_vode
		  }
	      }
	  }
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

    amrex::Print()<<"Maximum of repacked final solution: "<<S_old.max(5,0,0)<<
      "\nArkode max: "<<S_old.max(5,0,0)-2880041560931.35<<
      "\nArkode relative max: "<<(S_old.max(5,0,0)-2880041560931.35)/2880041560931.35<<
      "\nCvode max: "<<S_old.max(5,0,0)-2880345125756.03<<
      "\nCvode relative max: "<<(S_old.max(5,0,0)-2880345125756.03)/2880345125756.03<<std::endl;
    amrex::Print()<<"Minimum of repacked final solution: "<<S_old.min(5,0,0)<<
      "\nArkode min: "<<(S_old.min(5,0,0)-1829794600851.52)<<
      "\nArkode relative min: "<<(S_old.min(5,0,0)-1829794600851.52)/1829794600851.52<<
      "\nCvode min: "<<S_old.min(5,0,0)-1829795010808.2<<
      "\nCvode relative min: "<<(S_old.min(5,0,0)-1829795010808.2)/1829795010808.2<<std::endl;
    
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

#ifdef AMREX_USE_CUDA
template<class SparseLinearSolver, class SystemClass, size_t order>
__global__
void do_sdc_kernel(Real* y_initial, Real* y_final, 
		   Real start_time, Real end_time, Real start_timestep,
                   Real tolerance, size_t maximum_newton_iters, 
		   bool fail_if_maximum_newton, Real maximum_steps,
		   Real epsilon, size_t size, bool use_adaptive_timestep) {

  typedef SdcIntegrator<SparseLinearSolver,SystemClass,order> SdcIntClass;

  const size_t WarpBatchSize = 128;
  const size_t WarpSize = 32;
  size_t warp_batch_id = blockIdx.x * WarpBatchSize;
  size_t global_index, local_index;

  SystemClass ode_system;
  global_index = threadIdx.x + warp_batch_id;

  if (global_index >= size) return;

  for (local_index = threadIdx.x; local_index < WarpBatchSize && global_index < size; local_index += WarpSize) {
    global_index = local_index + warp_batch_id;

    SdcIntClass sdc;
    RealVector<SystemClass::neqs> y_ini;

    for (size_t i = 0; i < SystemClass::neqs; i++) {
      y_ini.data[i] = y_initial[global_index * SystemClass::neqs + i];
    }
    SdcIntClass::set_jacobian_layout(sdc, ode_system);
    printf("before evaluate y_fin=%g\n",y_ini.data[0]);
    SdcIntClass::initialize(sdc, y_ini, 
			    start_time, end_time, start_timestep,
			    tolerance, maximum_newton_iters, 
			    fail_if_maximum_newton, maximum_steps,
			    epsilon, use_adaptive_timestep);
    printf("before evaluate y_fin=%g\n",y_ini.data[0]);

    for (size_t i = 0; i < maximum_steps; i++) {
      SdcIntClass::prepare(sdc);
      printf("before evaluate y_fin=%g\n",y_ini.data[0]);
      SdcIntClass::solve(sdc);
      printf("before evaluate y_fin=%g\n",y_ini.data[0]);
      SdcIntClass::update(sdc);
      printf("before evaluate y_fin=%g\n",y_ini.data[0]);
      if (SdcIntClass::is_finished(sdc)) break;
    }
    printf("before evaluate y_fin=%g\n",y_ini.data[0]);
    RealVector<SystemClass::neqs>& y_fin = SdcIntClass::get_current_solution(sdc);
    printf("before evaluate y_fin=%g\n",y_fin.data[0]);

    for (size_t i = 0; i < SystemClass::neqs; i++) {
      y_final[global_index * SystemClass::neqs + i] = y_fin.data[i];
    }
  }
}
#endif

template<class SparseLinearSolver, class SystemClass, size_t order>
void do_sdc_host(Real* y_initial, Real* y_final, 
		 Real start_time, Real end_time, Real start_timestep,
		 Real tolerance, size_t maximum_newton_iters, 
		 bool fail_if_maximum_newton, Real maximum_steps,
		 Real epsilon, size_t size, bool use_adaptive_timestep) {

  typedef SdcIntegrator<SparseLinearSolver,SystemClass,order> SdcIntClass;

  SystemClass ode_system;

  for (size_t global_index = 0; global_index < size; global_index++) {
    SdcIntClass sdc;
    RealVector<SystemClass::neqs> y_ini;

    for (size_t i = 0; i < SystemClass::neqs; i++) {
      y_ini.data[i] = y_initial[global_index * SystemClass::neqs + i];
    }

    SdcIntClass::set_jacobian_layout(sdc, ode_system);
    SdcIntClass::initialize(sdc, y_ini, 
			    start_time, end_time, start_timestep,
			    tolerance, maximum_newton_iters, 
			    fail_if_maximum_newton, maximum_steps,
			    epsilon, use_adaptive_timestep);

    for (size_t i = 0; i < maximum_steps; i++) {
      SdcIntClass::prepare(sdc);
      SdcIntClass::solve(sdc);
      SdcIntClass::update(sdc);
      if (SdcIntClass::is_finished(sdc)) break;
    }

    RealVector<SystemClass::neqs>& y_fin = SdcIntClass::get_current_solution(sdc);
    for (size_t i = 0; i < SystemClass::neqs; i++) {
      y_final[global_index * SystemClass::neqs + i] = y_fin.data[i];
    }
  }
}


int integrate_state_vec_mfin
  (amrex::Array4<Real> const& state4,
   amrex::Array4<Real> const& diag_eos4,
   const Box& tbx,
   const Real& a, const Real& delta_time)
{

  realtype reltol, abstol;
  int flag;
    
  reltol = 1e-4;  /* Set the tolerances */
  abstol = 1e-4;

  int one_in = 1;
  

      const auto len = amrex::length(tbx);  // length of box
      const auto lo  = amrex::lbound(tbx);  // lower bound of box

      N_Vector u;
      void *cvode_mem;
      realtype t=0.0;
				
      u = NULL;
      cvode_mem = NULL;

      long int neq = len.x*len.y*len.z;
      amrex::Gpu::streamSynchronize();
      //      amrex::Print()<<"len: "<<len<<"lo: "<<lo<<"neq: "<<neq<<std::endl;
				int loop = 1;

#ifdef AMREX_USE_CUDA
			cudaStream_t currentStream = amrex::Gpu::Device::cudaStream();	
#ifndef MAKE_MANAGED

#ifdef PATCH
				u = N_VMakeWithManagedAllocator_Cuda(neq,sunalloc,sunfree);  /* Allocate u vector */
#else
				u = N_VNewManaged_Cuda(neq);  /* Allocate u vector */
#endif
				double* dptr=N_VGetDeviceArrayPointer_Cuda(u);

#ifdef PATCH
				N_Vector e_orig = N_VMakeWithManagedAllocator_Cuda(neq,sunalloc,sunfree);  /* Allocate u vector */
				N_Vector e_orig2 = N_VMakeWithManagedAllocator_Cuda(neq,sunalloc,sunfree);  /* Allocate u vector */
#else
				N_Vector e_orig = N_VNewManaged_Cuda(neq);  /* Allocate u vector */
				N_Vector e_orig2 = N_VNewManaged_Cuda(neq);  /* Allocate u vector */
#endif
				double* eptr=N_VGetDeviceArrayPointer_Cuda(e_orig);
				double* eptr2=N_VGetDeviceArrayPointer_Cuda(e_orig2);
				N_VSetCudaStream_Cuda(e_orig, &currentStream);
				N_VSetCudaStream_Cuda(u, &currentStream);

#ifdef PATCH
				N_Vector Data = N_VMakeWithManagedAllocator_Cuda(4*neq,sunalloc,sunfree);  // Allocate u vector 
#else
				N_Vector Data = N_VNewManaged_Cuda(4*neq);  /* Allocate u vector */
#endif
				double* rparh = N_VGetDeviceArrayPointer_Cuda(Data);
				N_VSetCudaStream_Cuda(Data, &currentStream);
				// shouldn't need to initialize 
				//N_VConst(0.0,Data);

#ifdef PATCH
				N_Vector abstol_vec = N_VMakeWithManagedAllocator_Cuda(neq,sunalloc,sunfree);  
#else
				N_Vector abstol_vec = N_VNewManaged_Cuda(neq);  /* Allocate u vector */
#endif
				double* abstol_ptr = N_VGetDeviceArrayPointer_Cuda(abstol_vec);
				N_VSetCudaStream_Cuda(abstol_vec,&currentStream);
				amrex::Gpu::Device::streamSynchronize();
#else
				double* dptr=(double*) The_Managed_Arena()->alloc(neq*sizeof(double));
				u = N_VMakeManaged_Cuda(neq,dptr);  /* Allocate u vector */
				double* eptr= (double*) The_Managed_Arena()->alloc(neq*sizeof(double));
				N_Vector e_orig = N_VMakeManaged_Cuda(neq,eptr);  /* Allocate u vector */
				N_VSetCudaStream_Cuda(e_orig, &currentStream);
				N_VSetCudaStream_Cuda(u, &currentStream);

				double* rparh = (double*) The_Managed_Arena()->alloc(4*neq*sizeof(double));
				N_Vector Data = N_VMakeManaged_Cuda(4*neq,rparh);  // Allocate u vector 
				N_VSetCudaStream_Cuda(Data, &currentStream);
				// shouldn't need to initialize 
				//N_VConst(0.0,Data);

				double* abstol_ptr = (double*) The_Managed_Arena()->alloc(neq*sizeof(double));
				N_Vector abstol_vec = N_VMakeManaged_Cuda(neq,abstol_ptr);
				N_VSetCudaStream_Cuda(abstol_vec,&currentStream);
#endif

#else
#ifdef _OPENMP
				int nthreads=omp_get_max_threads();
				u = N_VNew_OpenMP(neq,nthreads);  /* Allocate u vector */
				N_Vector e_orig = N_VNew_OpenMP(neq,nthreads);  /* Allocate u vector */
				double* eptr=N_VGetArrayPointer_Serial(e_orig);
				double* dptr=N_VGetArrayPointer_OpenMP(u);

				N_Vector Data = N_VNew_OpenMP(4*neq,nthreads);  // Allocate u vector 
				N_VConst(0.0,Data);
				double* rparh=N_VGetArrayPointer_OpenMP(Data);
				N_Vector abstol_vec = N_VNew_OpenMP(neq,nthreads);
				double* abstol_ptr=N_VGetArrayPointer_OpenMP(abstol_vec);
#else
				u = N_VNew_Serial(neq);  /* Allocate u vector */
				N_Vector e_orig = N_VNew_Serial(neq);  /* Allocate u vector */
				N_Vector e_orig2 = N_VNew_Serial(neq);  /* Allocate u vector */
				double* eptr=N_VGetArrayPointer_Serial(e_orig);
				double* eptr2=N_VGetArrayPointer_Serial(e_orig2);

				double* dptr=N_VGetArrayPointer_Serial(u);

				N_Vector Data = N_VNew_Serial(4*neq);  // Allocate u vector 
				N_VConst(0.0,Data);
				double* rparh=N_VGetArrayPointer_Serial(Data);
				N_Vector abstol_vec = N_VNew_Serial(neq);
				double* abstol_ptr=N_VGetArrayPointer_Serial(abstol_vec);
#endif
#endif

#ifdef _OPENMP
      const Dim3 hi = amrex::ubound(tbx);

#pragma omp parallel for collapse(3)
      for (int k = lo.z; k <= hi.z; ++k) {
	for (int j = lo.y; j <= hi.y; ++j) {
	    for (int i = lo.x; i <= hi.x; ++i) {
	      fort_ode_eos_setup(a,delta_time);
	      if(i==24&&j==14&&k==19)
		{
		  std::cout<<"rho_e(24,14,19): "<<state4(24,14,19,Eint)<<"\n"
			   <<state4(i,j,k,Eden)<<"\n"
			   <<state4(i,j,k,Density)<<"\n"
			   <<state4(i,j,k,Eint)/state4(i,j,k,Density)<<"\n"
			   <<diag_eos4(i,j,k,Ne_comp)<<"\n"
			   <<diag_eos4(i,j,k,Temp_comp)<<"\n"
			   <<"z: "<<1/a-1<<"\n"
			   <<"a: "<<a<<std::endl;
		}
#else
				AMREX_PARALLEL_FOR_3D ( tbx, i,j,k,
				{				  
#endif
				  int idx = i+j*len.x+k*len.x*len.y-(lo.x+lo.y*len.x+lo.z*len.x*len.y);
				  dptr[idx]=state4(i,j,k,Eint)/state4(i,j,k,Density);
				  eptr[idx]=state4(i,j,k,Eint)/state4(i,j,k,Density);
				  eptr2[idx]=state4(i,j,k,Eint)/state4(i,j,k,Density);
				  rparh[4*idx+0]= diag_eos4(i,j,k,Temp_comp);   //rpar(1)=T_vode
				  rparh[4*idx+1]= diag_eos4(i,j,k,Ne_comp);//    rpar(2)=ne_vode
				  rparh[4*idx+2]= state4(i,j,k,Density); //    rpar(3)=rho_vode
				  rparh[4*idx+3]=1/a-1;    //    rpar(4)=z_vode
				  //				  abstol_ptr[idx]= diag_eos4(i,j,k,Ne_comp)<1e-7||true ? state4(i,j,k,Eint)/state4(i,j,k,Density)*abstol : 1e4*state4(i,j,k,Eint)/state4(i,j,k,Density)*abstol ;
				  abstol_ptr[idx]= state4(i,j,k,6) ? state4(i,j,k,Eint)/state4(i,j,k,Density)*1e-12 : state4(i,j,k,Eint)/state4(i,j,k,Density)*abstol;
				  //				}
#ifdef _OPENMP
				}
				}
				}
#pragma omp barrier
#else
				});
      amrex::Gpu::Device::streamSynchronize();
#endif

      SUNMatrix A = NULL;             /* empty matrix for linear solver */
      SUNLinearSolver LS = NULL;      /* empty linear solver object */
     
      if(use_cvdiag==2)
	{
	  std::cout<<"Dense: "<<use_cvdiag<<std::endl;
	  A = SUNDenseMatrix(neq, neq);
	  LS = SUNLinSol_Dense(u, A);
	}
      else if(use_cvdiag==3||N_VMin(abstol_vec)<1e-7)
	{
	  LS = SUNLinSol_PCG(u,0,1);
	}
      else if(use_cvdiag==4)
	{
	  LS = SUNLinSol_SPGMR(u,0,1);
	  flag = SUNLinSol_SPGMRSetGSType(LS, CLASSICAL_GS);

	  //	  flag = ARKStepSetEpsLin(cvode_mem, DELT);
	  //	  if(check_flag(&flag, "ARKStepSetEpsLin", 1)) return(1);
	}
      amrex::Gpu::streamSynchronize();
  if(use_cvode==1)
  {
#ifdef CV_NEWTON
				cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
#else
				cvode_mem = CVodeCreate(CV_BDF);
#endif
				flag = CVodeInit(cvode_mem, f, t, u);

				//				N_VScale(abstol,u,abstol_vec);
				if(N_VMin(abstol_vec)<1e-7)
				  flag = CVodeSVtolerances(cvode_mem, 1e-6*reltol, abstol_vec);
				else
				  flag = CVodeSVtolerances(cvode_mem, reltol, abstol_vec);

				if(use_cvdiag==1)
				  flag = CVDiag(cvode_mem);
				else if(use_cvdiag==2)//||N_VMin(abstol_vec)<1e-7)
				  {
				  flag = CVodeSetLinearSolver(cvode_mem, LS, A);        /* Attach matrix and linear solver */
				  }
				else
				  {
				  flag = CVodeSetLinearSolver(cvode_mem, LS, NULL);        /* Attach matrix and linear solver */
				  }
				CVodeSetUserData(cvode_mem, &Data);
				CVodeSetMaxNumSteps(cvode_mem,2000);
				//				CVodeSetInitStep(cvode_mem, delta_time/5);
				flag = CVode(cvode_mem, delta_time, u, &t, CV_NORMAL);

  }

  else if(use_cvode==0)
    {
      
      cvode_mem = ARKStepCreate(NULL, f, t, u);

      //      N_VScale(abstol,u,abstol_vec);
      flag = ARKStepSVtolerances(cvode_mem, reltol, abstol_vec);
				//////				flag = CVodeSVtolerances(cvode_mem, reltol, abstol_vec);
				
      /* Linear solver interface */
      ARKStepSetLinear(cvode_mem,0);
      if(use_cvdiag==2)
	{
	  flag = ARKStepSetLinearSolver(cvode_mem, LS, A);        /* Attach matrix and linear solver */
	}
      else
	{
	  flag = ARKStepSetLinearSolver(cvode_mem, LS, NULL);        /* Attach matrix and linear solver */
	}
      flag = ARKStepSetAdaptivityMethod(cvode_mem, 4,1,0, NULL);
      //      flag = ARKStepSetInitStep(cvode_mem, delta_time/5);
		//				flag = CVodeSStolerances(cvode_mem, reltol, N_VMin(abstol_vec));
      //      flag = CVDiag(cvode_mem);

      ARKStepSetMaxNumSteps(cvode_mem,2000);
				/*
				N_Vector constrain=N_VClone(u);
				N_VConst(2,constrain);	      
				flag =CVodeSetConstraints(cvode_mem,constrain);
				*/
      ARKStepSetUserData(cvode_mem, &Data);
				//
				//				BL_PROFILE_VAR("Nyx::strang_second_cvode",cvode_timer2);
      flag = ARKStepEvolve(cvode_mem, delta_time, u, &t, ARK_NORMAL);
      amrex::Print()<<"tout"<<t<<std::endl;
				//				amrex::Gpu::Device::streamSynchronize();
				//				BL_PROFILE_VAR_STOP(cvode_timer2);

    }
  else
    {
  const size_t order = 4;
  Real start_time = 0.0;
  /*
  Real end_time = 1.6842553218299e-10;//delta_time;
  Real start_timestep = delta_time;(end_time - start_time)/5120.0;
  Real tolerance = 1.0e-1;*/
  Real end_time = delta_time;
  Real start_timestep = (end_time - start_time)/5120.0;
  start_timestep = delta_time;
  Real tolerance = 1.0e-12;
  size_t maximum_newton_iters = 1000;
  size_t maximum_steps = 1000000;
  bool fail_if_maximum_newton = true;
  Real epsilon = std::numeric_limits<Real>::epsilon();
  bool use_adaptive_timestep = true;
  size_t num_systems=neq/3;

  const int nThreads = 32;
  const size_t WarpBatchSize = 128;
  const int nBlocks = static_cast<int>(ceil(((double) num_systems)/(double) WarpBatchSize));

#ifndef AMREX_USE_CUDA
  do_sdc_host<SparseGaussJordan, VodeSystem, order>(eptr2, dptr,
						    start_time, end_time, start_timestep,
						    tolerance, maximum_newton_iters,
						    fail_if_maximum_newton, maximum_steps,
						    epsilon, num_systems, use_adaptive_timestep);
#else
  do_sdc_kernel<SparseGaussJordan, 
		VodeSystem, 
		order><<<nBlocks, nThreads>>>(eptr2, dptr,
					      start_time, end_time, start_timestep,
					      tolerance, maximum_newton_iters,
					      fail_if_maximum_newton, maximum_steps,
					      epsilon, num_systems, use_adaptive_timestep);

#endif

    }
      amrex::Gpu::streamSynchronize();
#ifdef AMREX_USE_CUDA
      N_VPrint_Cuda(u);
      N_VPrint_Cuda(e_orig2);
#else
      N_VPrint_Serial(u);
      N_VPrint_Serial(e_orig2);
#endif
#ifdef AMREX_DEBUG
				PrintFinalStats(cvode_mem);
#endif

				//				diag_eos(i,j,k,Temp_comp)=rparh[0];   //rpar(1)=T_vode
				//	diag_eos(i,j,k,Ne_comp)=rparh[1];//    rpar(2)=ne_vode
				// rho should not change  rho_tmp_ptr[i]=rparh[4*i+2]; //    rpar(3)=rho_vode
				/*
				const Dim3 hi = amrex::ubound(tbx);

				amrex::Gpu::LaunchSafeGuard lsg(true);
				amrex::Gpu::streamSynchronize();

				for (int k = lo.z; k <= hi.z; ++k) {
				  for (int j = lo.y; j <= hi.y; ++j) {
				    for (int i = lo.x; i <= hi.x; ++i) {
				      amrex::Gpu::streamSynchronize();
				      int idx = i+j*len.x+k*len.x*len.y-(lo.x+lo.y*len.x+lo.z*len.x*len.y);
				      realtype t=0.0;
				      amrex::Gpu::streamSynchronize();
				      if(rparh[4*idx+1] >1e-4)
					{
					  std::cout<<"rparh: "<<rparh[4*idx+1]<<std::endl;
					  N_Vector u = N_VNewManaged_Cuda(1);  
					  double* dptr_tmp=N_VGetDeviceArrayPointer_Cuda(u);
					  
					  N_Vector Data = N_VNewManaged_Cuda(4*1);  // Allocate u vector 
					  N_VConst(0.0,Data);

					  N_VSetCudaStream_Cuda(Data, &currentStream);
					  N_VSetCudaStream_Cuda(u, &currentStream);
					  double* rparh_tmp=N_VGetDeviceArrayPointer_Cuda(Data);
					  N_Vector abstol_vec = N_VNewManaged_Cuda(1);
					  N_VSetCudaStream_Cuda(abstol_vec, &currentStream);
					  double* abstol_ptr=N_VGetDeviceArrayPointer_Cuda(abstol_vec);
					  amrex::Gpu::streamSynchronize();
					  idx=0;

					  dptr_tmp[idx]=state4(i,j,k,Eint)/state4(i,j,k,Density);
					  rparh_tmp[4*idx+0]= diag_eos4(i,j,k,Temp_comp);   //rpar(1)=T_vode
					  rparh_tmp[4*idx+1]= diag_eos4(i,j,k,Ne_comp);//    rpar(2)=ne_vode
					  rparh_tmp[4*idx+2]= state4(i,j,k,Density); //    rpar(3)=rho_vode
					  rparh_tmp[4*idx+3]=1/a-1;    //    rpar(4)=z_vode
					  amrex::Gpu::streamSynchronize();
					  std::cout<<"solving new outputs"<<dptr_tmp[idx]<<"\t"<<rparh_tmp[4*idx+2]<<"\t"<<idx<<"\t"<<i<<"\t"<<j<<"\t"<<k<<std::endl;
#ifdef CV_NEWTON
					  void* cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
#else
					  void* cvode_mem = CVodeCreate(CV_BDF);
#endif
					  flag = CVodeInit(cvode_mem, f, t, u);

					  N_VScale(abstol,u,abstol_vec);

					  flag = CVodeSVtolerances(cvode_mem, reltol, abstol_vec);

				//				flag = CVodeSStolerances(cvode_mem, reltol, dptr[0]*abstol);
					  flag = CVDiag(cvode_mem);

					  CVodeSetMaxNumSteps(cvode_mem,2000);

					  CVodeSetUserData(cvode_mem, &Data);
				//				CVodeSetMaxStep(cvode_mem, delta_time/10);
				//				BL_PROFILE_VAR("Nyx::strang_second_cvode",cvode_timer2);
					  flag = CVode(cvode_mem, delta_time, u, &t, CV_NORMAL);
					  amrex::Gpu::streamSynchronize();

					  idx= i+j*len.x+k*len.x*len.y-(lo.x+lo.y*len.x+lo.z*len.x*len.y);
					  dptr[idx]=dptr_tmp[0];
					  rparh[4*idx+0]=rparh_tmp[0];
					  rparh[4*idx+1]=rparh_tmp[1];
					  rparh[4*idx+2]=rparh_tmp[2];
					  rparh[4*idx+3]=rparh_tmp[3];
					  std::cout<<"finished setting new outputs"<<dptr[idx]<<"\t"<<rparh[4*idx+2]<<idx<<i<<j<<k<<std::endl;

					}
				      amrex::Gpu::streamSynchronize();

				    }
				  }
				}

				std::cout<<"finished checking new outputs"<<std::endl;*/

#ifdef _OPENMP
#pragma omp parallel for collapse(3)
      for (int k = lo.z; k <= hi.z; ++k) {
	for (int j = lo.y; j <= hi.y; ++j) {
	    for (int i = lo.x; i <= hi.x; ++i) {
	      if(i==24&&j==14&&k==19)
		std::cout<<"rho_e(24,14,19): "<<state4(24,14,19,Eint)<<std::endl;
#else
				AMREX_PARALLEL_FOR_3D ( tbx, i,j,k,
				{				  
#endif
				  int  idx= i+j*len.x+k*len.x*len.y-(lo.x+lo.y*len.x+lo.z*len.x*len.y);
				//				for (int i= 0;i < neq; ++i) {
				  fort_ode_eos_finalize(&(dptr[idx*loop]), &(rparh[4*idx*loop]), one_in);
				  diag_eos4(i,j,k,Temp_comp)=rparh[4*idx*loop+0];   //rpar(1)=T_vode
				  diag_eos4(i,j,k,Ne_comp)=rparh[4*idx*loop+1];//    rpar(2)=ne_vode
				
				  state4(i,j,k,Eint)  += state4(i,j,k,Density) * (dptr[idx*loop]-eptr[idx]);
				  state4(i,j,k,Eden)  += state4(i,j,k,Density) * (dptr[idx*loop]-eptr[idx]);
				  //				}
				//PrintFinalStats(cvode_mem);
#ifdef _OPENMP
	      if(i==24&&j==14&&k==19)
		{
		std::cout<<"end: rho_e(24,14,19): "<<state4(24,14,19,Eint)<<std::endl;
		std::cout<<"end: Temp(24,14,19): "<<diag_eos4(24,14,19,Temp_comp)<<std::endl;
		std::cout<<"end: Ne(24,14,19): "<<diag_eos4(24,14,19,Ne_comp)<<std::endl;
		}
				}
				}
				}
#pragma omp barrier
#else
				});
      amrex::Gpu::Device::streamSynchronize();
#endif

#ifdef MAKE_MANAGED
#ifdef AMREX_USE_GPU
      The_Managed_Arena()->free(dptr);
      The_Managed_Arena()->free(eptr);
      //      The_Managed_Arena()->free(constrain);
      The_Managed_Arena()->free(rparh);
      The_Managed_Arena()->free(abstol_ptr);
#endif
#endif
				N_VDestroy(u);          /* Free the u vector */
				N_VDestroy(e_orig);          /* Free the e_orig vector */
				///			N_VDestroy(constrain);          /* Free the constrain vector */
				N_VDestroy(abstol_vec);          /* Free the u vector */
				N_VDestroy(Data);          /* Free the userdata vector */

				
				//N_Vector ele= N_VClone(u);
				if(use_cvode)
				  {
				    //				      CVodeGetEstLocalErrors(cvode_mem, ele);
				    //				      std::cout<<"maximum local estimated error: "<<N_VMaxNorm(ele)<<std::endl;
				    //				      amrex::Gpu::streamSynchronize();
				    //N_VDestroy(ele);          // Free the u vector 
				      CVodeFree(&cvode_mem);  // Free the integrator memory 
				  }
				else
				  {
				    //				    ARKStepGetEstLocalErrors(cvode_mem, ele);
				    //				    std::cout<<"maximum local estimated error: "<<N_VMaxNorm(ele)<<std::endl;
				    //				    amrex::Gpu::streamSynchronize();
				    //				    N_VDestroy(ele);          // Free the u vector 
				    ARKStepFree(&cvode_mem);  // Free the integrator memory 
				  }
			      //);
				/*			    }
			
	}
	}*/
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
  
      Real rpara[4];
      Real* rpar_ptr=&(rpara[0]);
      Real a = 0.0688707121;
      rpara[0]=11026.08482;
      rpara[1]=0.008728890037;
      rpara[2]=2.587879236e+10;
      rpara[3]=1/a-1;
  
  //  cudaStream_t currentStream = amrex::Cuda::Device::cudaStream();
#ifdef _OPENMP
  #pragma omp parallel for
  for(int idx=0;idx<neq;idx++)
      RhsFnReal(t,u_ptr+idx,udot_ptr+idx, rpar+4*idx, 1);
#else
  AMREX_LAUNCH_DEVICE_LAMBDA ( neq, idx, {
      RhsFnReal(t,u_ptr+idx,udot_ptr+idx, rpar_ptr, 1);
  });
#endif
  amrex::Gpu::Device::streamSynchronize();//cudaStreamSynchronize(currentStream);
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

/*-------------------------------
 * Private helper functions
 *-------------------------------*/

/* Check function return value...
    opt == 0 means SUNDIALS function allocates memory so check if
             returned NULL pointer
    opt == 1 means SUNDIALS function returns a flag so check if
             flag >= 0
    opt == 2 means function allocates memory so check if returned
             NULL pointer
*/
static int check_flag(void *flagvalue, const char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1; }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
              funcname, *errflag);
      return 1; }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1; }

  return 0;
}

static void PrintFinalStats(void *cvode_mem)
{
  long lenrw, leniw ;
  long lenrwLS, leniwLS;
  long int nst, nfe, nsetups, nni, ncfn, netf;
  long int nli, npe, nps, ncfl, nfeLS;
  long int nst_a, nfi, nje;
  int retval, flag;

  if(use_cvode==1)
    {
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
    }
  else if(use_cvode==0)
    {
  /* Get/print some final statistics on how the solve progressed */
  flag = ARKStepGetNumSteps(cvode_mem, &nst);
  check_flag(&flag, "ARKStepGetNumSteps", 1);
  flag = ARKStepGetNumStepAttempts(cvode_mem, &nst_a);
  check_flag(&flag, "ARKStepGetNumStepAttempts", 1);
  flag = ARKStepGetNumRhsEvals(cvode_mem, &nfe, &nfi);
  check_flag(&flag, "ARKStepGetNumRhsEvals", 1);
  flag = ARKStepGetNumLinSolvSetups(cvode_mem, &nsetups);
  check_flag(&flag, "ARKStepGetNumLinSolvSetups", 1);
  flag = ARKStepGetNumErrTestFails(cvode_mem, &netf);
  check_flag(&flag, "ARKStepGetNumErrTestFails", 1);
  flag = ARKStepGetNumNonlinSolvIters(cvode_mem, &nni);
  check_flag(&flag, "ARKStepGetNumNonlinSolvIters", 1);
  flag = ARKStepGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  check_flag(&flag, "ARKStepGetNumNonlinSolvConvFails", 1);
  flag = ARKStepGetNumJacEvals(cvode_mem, &nje);
  check_flag(&flag, "ARKStepGetNumJacEvals", 1);
  flag = ARKStepGetNumLinRhsEvals(cvode_mem, &nfeLS);
  check_flag(&flag, "ARKStepGetNumLinRhsEvals", 1);

  printf("\nFinal Solver Statistics:\n");
  printf("   Internal solver steps = %li (attempted = %li)\n", nst, nst_a);
  printf("   Total RHS evals:  Fe = %li,  Fi = %li\n", nfe, nfi);
  printf("   Total linear solver setups = %li\n", nsetups);
  printf("   Total RHS evals for setting up the linear system = %li\n", nfeLS);
  printf("   Total number of Jacobian evaluations = %li\n", nje);
  printf("   Total number of Newton iterations = %li\n", nni);
  printf("   Total number of linear solver convergence failures = %li\n", ncfn);
  printf("   Total number of error test failures = %li\n\n", netf);
    }

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
