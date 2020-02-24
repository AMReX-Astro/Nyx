#include <fstream>
#include <iomanip>

#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>
#include <AMReX_PlotFileUtil.H>

#include <AMReX_BLFort.H>
#include <Nyx.H>
#include <Nyx_F.H>

#include <cvode/cvode.h>               /* prototypes for CVODE fcts., consts. */
#include <cvode/cvode_diag.h>          /* access to CVDiag interface */
#include <sundials/sundials_types.h>   /* definition of type realtype */
#include <sundials/sundials_math.h>    /* definition of ABS and EXP   */

#include <nvector/nvector_serial.h>
#ifdef AMREX_USE_CUDA
#include <nvector/nvector_cuda.h>
#endif
using namespace amrex;

/* Functions Called by the Solver */
static int f(realtype t, N_Vector u, N_Vector udot, void *user_data);

static void PrintFinalStats(void *cvode_mem);

static void PrintOutput(realtype t, realtype umax, long int nst);

/* Private function to check function return values */
static int check_retval(void *flagvalue, const char *funcname, int opt);

int Nyx::integrate_state_cell
  (amrex::MultiFab &S_old,
   amrex::MultiFab &D_old,
   const Real& a, const Real& delta_time)
{
    // time = starting time in the simulation
  realtype reltol, abstol;
    
  reltol = 1e-4;  /* Set the tolerances */
  abstol = 1e-4;

  #ifdef _OPENMP
  #pragma omp parallel if (Gpu::notInLaunchRegion())
  #endif
  for ( MFIter mfi(S_old, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
      fort_ode_eos_setup(a,delta_time);
      //check that copy contructor vs create constructor works??
      const Box& tbx = mfi.tilebox();
      Array4<Real> const& state = S_old.array(mfi);
      Array4<Real> const& diag_eos = D_old.array(mfi);
      long int neq=1;
      const Dim3 lo = amrex::lbound(tbx);
      const Dim3 hi = amrex::ubound(tbx);

      N_Vector u;
      //  SUNLinearSolver LS;
      void *cvode_mem;
      realtype t=0.0;
      double* dptr;
      realtype tout;
      realtype umax;
      int flag;
      
      u = NULL;
      //  LS = NULL;
      cvode_mem = NULL;
      
      u = N_VNew_Serial(1);  /* Allocate u vector */
      dptr=N_VGetArrayPointer_Serial(u);
      //if(check_retval((void*)u, "N_VNew_Serial", 0)) return(1);
      
      N_Vector Data = N_VNew_Serial(4);  // Allocate u vector 
      N_VConst(0.0,Data);
      double* rparh=N_VGetArrayPointer_Serial(Data);
      
#ifdef CV_NEWTON
      cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
#else
      cvode_mem = CVodeCreate(CV_BDF);
#endif
      flag = CVodeInit(cvode_mem, f, t, u);
      flag = CVodeSStolerances(cvode_mem, reltol, dptr[0]*abstol);
      flag = CVDiag(cvode_mem);
      
      N_Vector constrain=N_VClone(u);
      N_VConst(2,constrain);	      
      flag =CVodeSetConstraints(cvode_mem,constrain);

      CVodeSetUserData(cvode_mem, &Data);
      for (int k = lo.z; k <= hi.z; ++k) {
	for (int j = lo.y; j <= hi.y; ++j) {
	  AMREX_PRAGMA_SIMD
	    for (int i = lo.x; i <= hi.x; ++i) {
	      t=0;
	      
	      N_VConst(state(i,j,k,Eint)/state(i,j,k,Density),u);
	      state(i,j,k,Eint)  -= state(i,j,k,Density) * (dptr[0]);
	      state(i,j,k,Eden)  -= state(i,j,k,Density) * (dptr[0]);
	      rparh[0]= diag_eos(i,j,k,Temp_comp);   //rpar(1)=T_vode
	      rparh[1]= diag_eos(i,j,k,Ne_comp);//    rpar(2)=ne_vode
	      rparh[2]= state(i,j,k,Density); //    rpar(3)=rho_vode
	      rparh[3]=1/a-1;    //    rpar(4)=z_vode
	      flag = CVodeReInit(cvode_mem, t, u);
	      flag = CVodeSStolerances(cvode_mem, reltol, dptr[0]*abstol);
	      
	      BL_PROFILE_VAR("Nyx::strang_second_cvode",cvode_timer2);
	      flag = CVode(cvode_mem, delta_time, u, &t, CV_NORMAL);
	      BL_PROFILE_VAR_STOP(cvode_timer2);
	      
	      AMREX_LAUNCH_DEVICE_LAMBDA(neq,i,
					 {
					   fort_ode_eos_finalize(&(dptr[0]), &(rparh[0]), 1);
					 });
	      diag_eos(i,j,k,Temp_comp)=rparh[0];   //rpar(1)=T_vode
	      diag_eos(i,j,k,Ne_comp)=rparh[1];//    rpar(2)=ne_vode
	      
	      state(i,j,k,Eint)  += state(i,j,k,Density) * (dptr[0]);
	      state(i,j,k,Eden)  += state(i,j,k,Density) * (dptr[0]);
	      
	    }
	}
      }
      N_VDestroy(u);          /* Free the u vector */
      N_VDestroy(Data);          /* Free the userdata vector */
      CVodeFree(&cvode_mem);  /* Free the integrator memory */

    }
  #ifdef AMREX_DEBUG
  #ifdef AMREX_DEBUG
        if (S_old.contains_nan())
            amrex::Abort("state has NaNs after the second strang call");
  #endif
#endif

    return 0;
}

int Nyx::integrate_state_growncell
  (amrex::MultiFab &S_old,
   amrex::MultiFab &D_old,
   const Real& a, const Real& delta_time)
{
    // time = starting time in the simulation
  realtype reltol, abstol;
  N_Vector u;
  //  SUNLinearSolver LS;
  int iout, flag;
  bool do_tiling=false;    

  reltol = 1e-4;  /* Set the tolerances */
  abstol = 1e-4;

  fort_ode_eos_setup(a,delta_time);
  amrex::Gpu::setLaunchRegion(false);
  #ifdef _OPENMP
  #pragma omp parallel if (Gpu::notInLaunchRegion())
  #endif
  for ( MFIter mfi(S_old, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
  fort_ode_eos_setup(a,delta_time);
      //check that copy contructor vs create constructor works??
      const Box& tbx = mfi.growntilebox();
      Array4<Real> const& state = S_old.array(mfi);
      Array4<Real> const& diag_eos = D_old.array(mfi);

      double rho,e_orig;
      long int neq=1;
      const Dim3 lo = amrex::lbound(tbx);
      const Dim3 hi = amrex::ubound(tbx);

      N_Vector u;
      //  SUNLinearSolver LS;
      void *cvode_mem;
      realtype t=0.0;
      double* dptr;
      realtype tout;
      realtype umax;
      int flag;
				
      u = NULL;
      //  LS = NULL;
      cvode_mem = NULL;
      
      u = N_VNew_Serial(1);  /* Allocate u vector */
      dptr=N_VGetArrayPointer_Serial(u);
      //if(check_retval((void*)u, "N_VNew_Serial", 0)) return(1);
      
      N_Vector Data = N_VNew_Serial(4);  // Allocate u vector 
      N_VConst(0.0,Data);
      double* rparh=N_VGetArrayPointer_Serial(Data);
      
#ifdef CV_NEWTON
      cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
#else
      cvode_mem = CVodeCreate(CV_BDF);
#endif
      flag = CVodeInit(cvode_mem, f, t, u);
      flag = CVodeSStolerances(cvode_mem, reltol, dptr[0]*abstol);
      flag = CVDiag(cvode_mem);

      N_Vector constrain=N_VClone(u);
      N_VConst(2,constrain);	      
      flag =CVodeSetConstraints(cvode_mem,constrain);
      
      CVodeSetUserData(cvode_mem, &Data);
      for (int k = lo.z; k <= hi.z; ++k) {
	for (int j = lo.y; j <= hi.y; ++j) {
	  AMREX_PRAGMA_SIMD
	    for (int i = lo.x; i <= hi.x; ++i) {
	      t=0;
	      
	      N_VConst(state(i,j,k,Eint)/state(i,j,k,Density),u);
	      state(i,j,k,Eint)  -= state(i,j,k,Density) * (dptr[0]);
	      state(i,j,k,Eden)  -= state(i,j,k,Density) * (dptr[0]);
	      rparh[0]= diag_eos(i,j,k,Temp_comp);   //rpar(1)=T_vode
	      rparh[1]= diag_eos(i,j,k,Ne_comp);//    rpar(2)=ne_vode
	      rparh[2]= state(i,j,k,Density); //    rpar(3)=rho_vode
	      rparh[3]=1/a-1;    //    rpar(4)=z_vode
	      flag = CVodeReInit(cvode_mem, t, u);
	      flag = CVodeSStolerances(cvode_mem, reltol, dptr[0]*abstol);
	      
	      BL_PROFILE_VAR("Nyx::strang_first_cvode",cvode_timer1);
	      flag = CVode(cvode_mem, delta_time, u, &t, CV_NORMAL);
	      BL_PROFILE_VAR_STOP(cvode_timer1);
	      
	      AMREX_LAUNCH_DEVICE_LAMBDA(neq,i,
					 {
					   fort_ode_eos_finalize(&(dptr[0]), &(rparh[0]), 1);
					 });
	      diag_eos(i,j,k,Temp_comp)=rparh[0];   //rpar(1)=T_vode
	      diag_eos(i,j,k,Ne_comp)=rparh[1];//    rpar(2)=ne_vode
	      
	      state(i,j,k,Eint)  += state(i,j,k,Density) * (dptr[0]);
	      state(i,j,k,Eden)  += state(i,j,k,Density) * (dptr[0]);
	      
	    }
	}
      }
      N_VDestroy(u);          /* Free the u vector */
      N_VDestroy(Data);          /* Free the userdata vector */
      CVodeFree(&cvode_mem);  /* Free the integrator memory */
      
    }
  #ifdef AMREX_DEBUG
        if (S_old.contains_nan())
            amrex::Abort("state has NaNs after the first strang call");
  #endif

    return 0;
}

#ifdef AMREX_USE_CUDA
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
  cudaStream_t currentStream = amrex::Gpu::Device::cudaStream();
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
  AMREX_LAUNCH_DEVICE_LAMBDA ( neq, idx, {
      //  f_rhs_test(t,u_ptr,udot_ptr, rpar, neq);
      RhsFnReal(t,u_ptr+idx,udot_ptr+idx, rpar+4*idx, 1);
  });
  cudaStreamSynchronize(currentStream);
  AMREX_GPU_ERROR_CHECK();

/*      N_VCopyFromDevice_Cuda(*(static_cast<N_Vector*>(user_data)));
      fprintf(stdout,"\nafter rparh[0]=%g \n\n",rparh[0]);
  fprintf(stdout,"\nafter rparh[1]=%g \n\n",rparh[1]);
  fprintf(stdout,"\nafter rparh[2]=%g \n\n",rparh[2]);
  fprintf(stdout,"\nafter rparh[3]=%g \n\n",rparh[3]);
  fprintf(stdout,"\nafter last rparh[4*(neq-1)+1]=%g \n\n",rparh[4*(neq-1)+1]);*/
  return 0;
}

#else
static int f(realtype t, N_Vector u, N_Vector udot, void* user_data)
{

  Real* udot_ptr=N_VGetArrayPointer_Serial(udot);
  Real* u_ptr=N_VGetArrayPointer_Serial(u);
  int neq=N_VGetLength_Serial(udot);
  double*  rpar=N_VGetArrayPointer_Serial(*(static_cast<N_Vector*>(user_data)));
  /*   fprintf(stdout,"\nt=%g \n\n",t);
  fprintf(stdout,"\nrparh[0]=%g \n\n",rpar[0]);
  fprintf(stdout,"\nrparh[1]=%g \n\n",rpar[1]);
  fprintf(stdout,"\nrparh[2]=%g \n\n",rpar[2]);
  fprintf(stdout,"\nrparh[3]=%g \n\n",rpar[3]);*/
  for(int tid=0;tid<neq;tid++)
    {
      //    fprintf(stdout,"\nrpar[4*tid+0]=%g\n",rpar[4*tid]);
    RhsFnReal(t,&(u_ptr[tid]),&(udot_ptr[tid]),&(rpar[4*tid]),1);
    //    fprintf(stdout,"\nafter rpar[4*tid+0]=%g\n",rpar[4*tid]);
    }
  /*      fprintf(stdout,"\nafter rparh[0]=%g \n\n",rpar[0]);
  fprintf(stdout,"\nafter rparh[1]=%g \n\n",rpar[1]);
  fprintf(stdout,"\nafter rparh[2]=%g \n\n",rpar[2]);
  fprintf(stdout,"\nafter rparh[3]=%g \n\n",rpar[3]);
  fprintf(stdout,"\nafter last rparh[4*(neq-1)+1]=%g \n\n",rpar[4*(neq-1)+1]);*/
  return 0;
}
#endif


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

  if (ParallelDescriptor::IOProcessor())
    {
  printf("\nFinal Statistics.. \n\n");
  printf("lenrw   = %5ld     leniw   = %5ld\n"  , lenrw, leniw);
  printf("nst     = %5ld\n"                     , nst);
  printf("nfe     = %5ld     nfeLS   = %5ld\n"  , nfe, nfeLS);
  printf("nsetups = %5ld     netf    = %5ld\n"  , nsetups, netf);
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
