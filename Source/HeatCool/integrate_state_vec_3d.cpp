#include <AMReX_Config.H>
#include <fstream>
#include <iomanip>

#include <cvode/cvode.h>               /* prototypes for CVODE fcts., consts. */
#include <cvode/cvode_diag.h>          /* access to CVDiag interface */
#include <sundials/sundials_types.h>   /* definition of type realtype */

#include <nvector/nvector_serial.h>
#ifndef AMREX_USE_GPU
#ifdef _OPENMP
#include <nvector/nvector_openmp.h>
#endif
#endif
#ifdef AMREX_USE_CUDA
#include <nvector/nvector_cuda.h>
#endif

#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Sundials.H>

#include <AMReX_BLFort.H>
#include <Nyx.H>
#include <f_rhs.H>

#define PATCH 1
//#define MAKE_MANAGED 1
using namespace amrex;

/* Functions Called by the Solver */
static int f(realtype t, N_Vector u, N_Vector udot, void *user_data);

static void PrintFinalStats(void *cvode_mem);

static void PrintOutput(realtype t, realtype umax, long int nst);

/* Private function to check function return values */
static int check_retval(void *flagvalue, const char *funcname, int opt);

int Nyx::integrate_state_vec
  (amrex::MultiFab &S_old,
   amrex::MultiFab &D_old,
   const Real& a, const Real& delta_time)
{
    // time = starting time in the simulation

  //Skip setup since parameters are hard-coded
  //  fort_ode_eos_setup(a,delta_time);
  long int store_steps=new_max_sundials_steps;

  bool tiling = (sundials_use_tiling && TilingIfNotGPU());
  for ( MFIter mfi(S_old, tiling ? MFItInfo().EnableTiling(sundials_tile_size) : MFItInfo()); mfi.isValid(); ++mfi )
  {

      //check that copy contructor vs create constructor works??
      const Box& tbx = mfi.tilebox();

      S_old[mfi].prefetchToDevice();
      D_old[mfi].prefetchToDevice();
      Array4<Real> const& state4 = S_old.array(mfi);
      Array4<Real> const& diag_eos4 = D_old.array(mfi);

      integrate_state_vec_mfin(state4,diag_eos4,tbx,a,delta_time,store_steps,new_max_sundials_steps);
  }
      return 0;
}

int Nyx::integrate_state_vec_mfin
  (amrex::Array4<Real> const& state4,
   amrex::Array4<Real> const& diag_eos4,
   const Box& tbx,
   const Real& a, const Real& delta_time,
   long int& old_max_steps, long int& new_max_steps)
{

  realtype reltol, abstol;
  int flag;

  reltol = sundials_reltol;  /* Set the tolerances */
  abstol = sundials_abstol;

  int one_in = 1;


      const auto len = amrex::length(tbx);  // length of box
      const auto lo  = amrex::lbound(tbx);  // lower bound of box

      N_Vector u;
      N_Vector e_orig;
      N_Vector Data;
      N_Vector abstol_vec;

      void *cvode_mem;
      realtype *dptr, *eptr, *rpar, *rparh, *abstol_ptr;
      realtype t=0.0;

      u = NULL;
      e_orig = NULL;
      Data = NULL;
      abstol_vec = NULL;
      cvode_mem = NULL;

      long int neq = len.x*len.y*len.z;
      amrex::Gpu::streamSynchronize();
      int loop = 1;

#ifdef AMREX_USE_CUDA
                        cudaStream_t currentStream = amrex::Gpu::Device::cudaStream();
                        //Choosing 256 here since this mimics Sundials default
                        // Possibly attempt to match n_threads_and_blocks
                        //     long N = ((tbx.numPts()+AMREX_GPU_NCELLS_PER_THREAD-1)/AMREX_GPU_NCELLS_PER_THREAD);
                        //     SUNCudaThreadDirectExecPolicy stream_exec_policy(AMREX_GPU_MAX_THREADS, currentStream);
                        //     SUNCudaGridStrideExecPolicy grid_exec_policy(AMREX_GPU_MAX_THREADS, std::max((N + AMREX_GPU_MAX_THREADS - 1) / AMREX_GPU_MAX_THREADS, static_cast<Long>(1)), currentStream);
                        //     SUNCudaBlockReduceExecPolicy reduce_exec_policy(AMREX_GPU_MAX_THREADS, std::max((N + AMREX_GPU_MAX_THREADS - 1) / AMREX_GPU_MAX_THREADS, static_cast<Long>(1)), currentStream);

                        SUNCudaThreadDirectExecPolicy stream_exec_policy(256, currentStream);
                        SUNCudaBlockReduceExecPolicy reduce_exec_policy(256, 0, currentStream);
                        if(sundials_alloc_type%2==0)
                        {
                          if(sundials_alloc_type==0)
                            u = N_VNewWithMemHelp_Cuda(neq, 1, *amrex::sundials::The_SUNMemory_Helper(),
                                  *amrex::sundials::The_Sundials_Context());  /* Allocate u vector */
                          else
                            u = N_VNewManaged_Cuda(neq, *amrex::sundials::The_Sundials_Context());  /* Allocate u vector */

                          dptr=N_VGetDeviceArrayPointer_Cuda(u);

                          if(sundials_alloc_type==0)
                            e_orig = N_VNewWithMemHelp_Cuda(neq, 1, *amrex::sundials::The_SUNMemory_Helper(),
                                      *amrex::sundials::The_Sundials_Context());  /* Allocate u vector */
                          else
                            e_orig = N_VNewManaged_Cuda(neq, *amrex::sundials::The_Sundials_Context());  /* Allocate u vector */

                          eptr=N_VGetDeviceArrayPointer_Cuda(e_orig);
                          N_VSetKernelExecPolicy_Cuda(e_orig, &stream_exec_policy, &reduce_exec_policy);
                          N_VSetKernelExecPolicy_Cuda(u, &stream_exec_policy, &reduce_exec_policy);

                          if(sundials_alloc_type==0)
                            Data = N_VNewWithMemHelp_Cuda(4*neq, 1, *amrex::sundials::The_SUNMemory_Helper(),
                                    *amrex::sundials::The_Sundials_Context());  // Allocate u vector
                          else
                            Data = N_VNewManaged_Cuda(4*neq, *amrex::sundials::The_Sundials_Context());  /* Allocate u vector */

                          rparh = N_VGetDeviceArrayPointer_Cuda(Data);
                          N_VSetKernelExecPolicy_Cuda(Data, &stream_exec_policy, &reduce_exec_policy);
                          // shouldn't need to initialize
                          //N_VConst(0.0,Data);

                          if(sundials_alloc_type==0)
                            abstol_vec = N_VNewWithMemHelp_Cuda(neq, 1, *amrex::sundials::The_SUNMemory_Helper(),
                                          *amrex::sundials::The_Sundials_Context());
                          else
                            abstol_vec = N_VNewManaged_Cuda(neq, *amrex::sundials::The_Sundials_Context());  /* Allocate u vector */

                          abstol_ptr = N_VGetDeviceArrayPointer_Cuda(abstol_vec);
                          N_VSetKernelExecPolicy_Cuda(abstol_vec, &stream_exec_policy, &reduce_exec_policy);
                          amrex::Gpu::Device::streamSynchronize();
                        }
                        else
                        {
                          dptr=(realtype*) The_Arena()->alloc(neq*sizeof(realtype));
                          u = N_VMakeManaged_Cuda(neq,dptr, *amrex::sundials::The_Sundials_Context());  /* Allocate u vector */
                          eptr= (realtype*) The_Arena()->alloc(neq*sizeof(realtype));
                          e_orig = N_VMakeManaged_Cuda(neq,eptr, *amrex::sundials::The_Sundials_Context());  /* Allocate u vector */
                          N_VSetKernelExecPolicy_Cuda(e_orig, &stream_exec_policy, &reduce_exec_policy);
                          N_VSetKernelExecPolicy_Cuda(u, &stream_exec_policy, &reduce_exec_policy);

                          rparh = (realtype*) The_Arena()->alloc(4*neq*sizeof(realtype));
                          Data = N_VMakeManaged_Cuda(4*neq,rparh, *amrex::sundials::The_Sundials_Context());  // Allocate u vector
                          N_VSetKernelExecPolicy_Cuda(Data, &stream_exec_policy, &reduce_exec_policy);
                          // shouldn't need to initialize
                          //N_VConst(0.0,Data);

                          abstol_ptr = (realtype*) The_Arena()->alloc(neq*sizeof(realtype));
                          abstol_vec = N_VMakeManaged_Cuda(neq,abstol_ptr, *amrex::sundials::The_Sundials_Context());
                          N_VSetKernelExecPolicy_Cuda(abstol_vec, &stream_exec_policy, &reduce_exec_policy);
                          amrex::Gpu::streamSynchronize();
                        }

#else
#if (defined(_OPENMP) && !defined(AMREX_USE_GPU))
                        int nthreads=omp_get_max_threads();
                        u = N_VNew_OpenMP(neq,nthreads, *amrex::sundials::The_Sundials_Context());  /* Allocate u vector */
                        e_orig = N_VNew_OpenMP(neq,nthreads, *amrex::sundials::The_Sundials_Context());  /* Allocate u vector */
                        eptr=N_VGetArrayPointer_Serial(e_orig);
                        dptr=N_VGetArrayPointer_OpenMP(u);

                        Data = N_VNew_OpenMP(4*neq,nthreads, *amrex::sundials::The_Sundials_Context());  // Allocate u vector
                        N_VConst(0.0,Data);
                        rparh=N_VGetArrayPointer_OpenMP(Data);
                        abstol_vec = N_VNew_OpenMP(neq,nthreads, *amrex::sundials::The_Sundials_Context());
                        abstol_ptr=N_VGetArrayPointer_OpenMP(abstol_vec);
#else
                        u = N_VNew_Serial(neq, *amrex::sundials::The_Sundials_Context());  /* Allocate u vector */
                        e_orig = N_VNew_Serial(neq, *amrex::sundials::The_Sundials_Context());  /* Allocate u vector */
                        eptr=N_VGetArrayPointer_Serial(e_orig);
                        dptr=N_VGetArrayPointer_Serial(u);

                        Data = N_VNew_Serial(4*neq, *amrex::sundials::The_Sundials_Context());  // Allocate u vector
                        N_VConst(0.0,Data);
                        rparh=N_VGetArrayPointer_Serial(Data);
                        abstol_vec = N_VNew_Serial(neq, *amrex::sundials::The_Sundials_Context());
                        abstol_ptr=N_VGetArrayPointer_Serial(abstol_vec);
#endif
#endif

#ifdef AMREX_USE_GPU
                                AMREX_PARALLEL_FOR_3D ( tbx, i,j,k,
                                {
#else
#ifdef _OPENMP
      const Dim3 hi = amrex::ubound(tbx);
#pragma omp parallel for collapse(3)
      for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                //Skip setup since parameters are hard-coded
#else
                                AMREX_PARALLEL_FOR_3D ( tbx, i,j,k,
                                {
#endif
#endif
                                  int idx = i+j*len.x+k*len.x*len.y-(lo.x+lo.y*len.x+lo.z*len.x*len.y);
                                  dptr[idx]=state4(i,j,k,Eint_comp)/state4(i,j,k,Density_comp);
                                  eptr[idx]=state4(i,j,k,Eint_comp)/state4(i,j,k,Density_comp);
                                  rparh[4*idx+0]= diag_eos4(i,j,k,Temp_comp);   //rpar(1)=T_vode
                                  rparh[4*idx+1]= diag_eos4(i,j,k,Ne_comp);//    rpar(2)=ne_vode
                                  rparh[4*idx+2]= state4(i,j,k,Density_comp); //    rpar(3)=rho_vode
                                  rparh[4*idx+3]=1/a-1;    //    rpar(4)=z_vode
                                  abstol_ptr[idx]= state4(i,j,k,Eint_comp)/state4(i,j,k,Density_comp)*abstol;
                                  //                            }
#ifdef AMREX_USE_GPU
                                });
                                amrex::Gpu::Device::streamSynchronize();
#else
#ifdef _OPENMP
                                }
                                }
                                }
#pragma omp barrier
#else
                                });
                                amrex::Gpu::Device::streamSynchronize();
#endif
#endif

                                cvode_mem = CVodeCreate(CV_BDF, *amrex::sundials::The_Sundials_Context());
                                flag = CVodeInit(cvode_mem, f, t, u);

                                N_VScale(abstol,u,abstol_vec);
                                //                              N_VConst(N_VMin(abstol_vec),abstol_vec);

                                flag = CVodeSVtolerances(cvode_mem, reltol, abstol_vec);

                                //                              flag = CVodeSStolerances(cvode_mem, reltol, dptr[0]*abstol);
                                flag = CVDiag(cvode_mem);

                                CVodeSetMaxNumSteps(cvode_mem,2000);

                                if(use_typical_steps)
                                    CVodeSetMaxStep(cvode_mem,delta_time/(old_max_steps));

                                N_Vector constrain;
                                if(use_sundials_constraint)
                                  {
                                    constrain=N_VClone(u);
                                    N_VConst(2,constrain);
                                    flag =CVodeSetConstraints(cvode_mem,constrain);
                                  }

#ifdef SUNDIALS_BUILD_PACKAGE_FUSED_KERNELS
                                if(use_sundials_fused)
                                {
                                     flag = CVodeSetUseIntegratorFusedKernels(cvode_mem, SUNTRUE);
                                }
#endif
                                CVodeSetUserData(cvode_mem, &Data);
                                //                              CVodeSetMaxStep(cvode_mem, delta_time/10);
                                //                              BL_PROFILE_VAR("Nyx::strang_second_cvode",cvode_timer2);
                                flag = CVode(cvode_mem, delta_time, u, &t, CV_NORMAL);
                                if(use_typical_steps)
                                  {
                                    long int nst=0;
                                    flag = CVodeGetNumSteps(cvode_mem, &nst);
                                    new_max_steps=std::max(nst,new_max_steps);
                                  }
                                //                              amrex::Gpu::Device::streamSynchronize();
                                //                              BL_PROFILE_VAR_STOP(cvode_timer2);

#ifdef AMREX_DEBUG
                                PrintFinalStats(cvode_mem);
#endif
  auto atomic_rates = atomic_rates_glob;
  Real lh_species = Nyx::h_species;
#ifdef AMREX_USE_GPU
                                AMREX_PARALLEL_FOR_3D ( tbx, i,j,k,
                                {
#else
#ifdef _OPENMP
#pragma omp parallel for collapse(3)
      for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                //Skip setup since parameters are hard-coded
#else
                                AMREX_PARALLEL_FOR_3D ( tbx, i,j,k,
                                {
#endif
#endif
                                  int  idx= i+j*len.x+k*len.x*len.y-(lo.x+lo.y*len.x+lo.z*len.x*len.y);
                                //                              for (int i= 0;i < neq; ++i) {
                                  ode_eos_finalize((dptr[idx*loop]), &(rparh[4*idx*loop]), one_in, atomic_rates, lh_species);
                                  diag_eos4(i,j,k,Temp_comp)=rparh[4*idx*loop+0];   //rpar(1)=T_vode
                                  diag_eos4(i,j,k,Ne_comp)=rparh[4*idx*loop+1];//    rpar(2)=ne_vode

                                  state4(i,j,k,Eint_comp)  += state4(i,j,k,Density_comp) * (dptr[idx*loop]-eptr[idx]);
                                  state4(i,j,k,Eden_comp)  += state4(i,j,k,Density_comp) * (dptr[idx*loop]-eptr[idx]);
                                  //                            }
                                //PrintFinalStats(cvode_mem);

#ifdef AMREX_USE_GPU
                                });
                                amrex::Gpu::Device::streamSynchronize();
#else
#ifdef _OPENMP
                                }
                                }
                                }
#pragma omp barrier
#else
                                });
                                amrex::Gpu::Device::streamSynchronize();
#endif
#endif

#ifdef AMREX_USE_CUDA
      if(sundials_alloc_type%2!=0)
      {
        The_Arena()->free(dptr);
        The_Arena()->free(eptr);
        if(use_sundials_constraint)
          The_Arena()->free(constrain);
        The_Arena()->free(rparh);
        The_Arena()->free(abstol_ptr);
      }
#endif

                                N_VDestroy(u);          /* Free the u vector */
                                N_VDestroy(e_orig);          /* Free the e_orig vector */
                                if(use_sundials_constraint)
                                  N_VDestroy(constrain);          /* Free the constrain vector */
                                N_VDestroy(abstol_vec);          /* Free the u vector */
                                N_VDestroy(Data);          /* Free the userdata vector */
                                CVodeFree(&cvode_mem);  /* Free the integrator memory */
                              //);
                                /*                          }

        }
        }*/
                                return 0;
}

int Nyx::integrate_state_grownvec
  (amrex::MultiFab &S_old,
   amrex::MultiFab &D_old,
   const Real& a, const Real& delta_time)
{
    //Skip setup since parameters are hard-coded
    //  fort_ode_eos_setup(a,delta_time);


    long int store_steps=old_max_sundials_steps;

    const Real prev_time     = state[State_Type].prevTime();

    bool tiling = (sundials_use_tiling && TilingIfNotGPU());
    for ( MFIter mfi(S_old, tiling); mfi.isValid(); ++mfi )
    {
        //check that copy contructor vs create constructor works??
        const Box& tbx = mfi.growntilebox();

        S_old[mfi].prefetchToDevice();
        D_old[mfi].prefetchToDevice();

        Array4<Real> const& state4 = S_old.array(mfi);
        Array4<Real> const& diag_eos4 = D_old.array(mfi);

        integrate_state_vec_mfin(state4,diag_eos4,tbx,a,delta_time,store_steps,old_max_sundials_steps);
    }

    return 0;
}

#ifdef AMREX_USE_GPU
static int f(realtype t, N_Vector u, N_Vector udot, void *user_data)
{
  Real* udot_ptr=N_VGetDeviceArrayPointer(udot);
  Real* u_ptr=N_VGetDeviceArrayPointer(u);
  int neq=N_VGetLength(udot);
  realtype*  rpar=N_VGetDeviceArrayPointer(*(static_cast<N_Vector*>(user_data)));

  auto atomic_rates = atomic_rates_glob;
  Real lh_species = Nyx::h_species;
  AMREX_LAUNCH_DEVICE_LAMBDA ( neq, idx, {
      //  f_rhs_test(t,u_ptr,udot_ptr, rpar, neq);
      f_rhs_rpar(t,*(u_ptr+idx),*(udot_ptr+idx), (rpar+4*idx), atomic_rates, lh_species);
  });
  amrex::Gpu::streamSynchronize();

  return 0;
}

#else
static int f(realtype t, N_Vector u, N_Vector udot, void* user_data)
{

  Real* udot_ptr=N_VGetArrayPointer_Serial(udot);
  Real* u_ptr=N_VGetArrayPointer_Serial(u);
  int neq=N_VGetLength_Serial(udot);
  realtype*  rpar=N_VGetArrayPointer_Serial(*(static_cast<N_Vector*>(user_data)));
  auto atomic_rates = atomic_rates_glob;
  Real lh_species = Nyx::h_species;
  #pragma omp parallel for
  for(int tid=0;tid<neq;tid++)
    {
        f_rhs_rpar(t, (u_ptr[tid]),(udot_ptr[tid]),&(rpar[4*tid]), atomic_rates, lh_species);
    }

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

