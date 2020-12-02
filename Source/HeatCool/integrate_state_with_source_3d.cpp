#include <fstream>
#include <iomanip>

#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>
#include <AMReX_PlotFileUtil.H>

#include <AMReX_BLFort.H>
#include <Nyx.H>
#include <f_rhs.H>
#include <f_rhs_struct.H>

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
//#define MAKE_MANAGED 1
using namespace amrex;

/* Functions Called by the Solver */
static int f(realtype t, N_Vector u, N_Vector udot, void *user_data);

int Nyx::integrate_state_struct
  (amrex::MultiFab &S_old,
   amrex::MultiFab &S_new,
   amrex::MultiFab &D_old,
   amrex::MultiFab &hydro_src,
   amrex::MultiFab &IR,
   amrex::MultiFab &reset_src,
   const Real& a, const amrex::Real& a_end,
   const Real& delta_time,
   const int sdc_iter)
{
    // time = starting time in the simulation

  amrex::Gpu::LaunchSafeGuard lsg(true);
  long int store_steps=new_max_sundials_steps;

  //#ifdef _OPENMP
  //#pragma omp parallel if (Gpu::notInLaunchRegion())
  //#endif
#ifdef _OPENMP
  for ( MFIter mfi(S_old, false); mfi.isValid(); ++mfi )
#else
#ifdef AMREX_USE_GPU
    for ( MFIter mfi(S_old, MFItInfo()); mfi.isValid(); ++mfi)
#else
  for ( MFIter mfi(S_old, true); mfi.isValid(); ++mfi )
#endif
#endif
    {

      //check that copy contructor vs create constructor works??
      const Box& tbx = mfi.tilebox();

      S_old[mfi].prefetchToDevice();
      D_old[mfi].prefetchToDevice();
      S_new[mfi].prefetchToDevice();
      hydro_src[mfi].prefetchToDevice();
      reset_src[mfi].prefetchToDevice();
      IR[mfi].prefetchToDevice();
      Array4<Real> const& state4 = S_old.array(mfi);
      Array4<Real> const& diag_eos4 = D_old.array(mfi);
      Array4<Real> const& state_n4 = S_new.array(mfi);
      Array4<Real> const& hydro_src4 = hydro_src.array(mfi);
      Array4<Real> const& reset_src4 = reset_src.array(mfi);
      Array4<Real> const& IR4 = IR.array(mfi);

      integrate_state_struct_mfin(state4,diag_eos4,state_n4,hydro_src4,reset_src4,IR4,tbx,a,a_end,delta_time,store_steps,new_max_sundials_steps,sdc_iter);
    }
    return 0;
}

int Nyx::integrate_state_struct_mfin
  (amrex::Array4<Real> const& state4,
   amrex::Array4<Real> const& diag_eos4,
   amrex::Array4<Real> const& state_n4,
   amrex::Array4<Real> const& hydro_src4,
   amrex::Array4<Real> const& reset_src4,
   amrex::Array4<Real> const& IR4,
   const Box& tbx,
   const Real& a, const amrex::Real& a_end, const Real& delta_time,
   long int& old_max_steps, long int& new_max_steps,
   const int sdc_iter)
{

    auto atomic_rates = atomic_rates_glob;
    auto f_rhs_data = (RhsData*)The_Arena()->alloc(sizeof(RhsData));
    Real gamma_minus_1 = gamma - 1.0;
    ode_eos_setup(f_rhs_data, gamma_minus_1, h_species);

    realtype reltol, abstol;
    int flag;
    
    reltol = 1e-4;  /* Set the tolerances */
    abstol = 1e-4;

    int one_in = 1;

      const auto len = amrex::length(tbx);  // length of box
      const auto lo  = amrex::lbound(tbx);  // lower bound of box

      N_Vector u;
      N_Vector e_orig;
      N_Vector Data;
      N_Vector abstol_vec;
      N_Vector T_vec;
      N_Vector ne_vec;
      N_Vector rho_vec;
      N_Vector rho_init_vec;
      N_Vector rho_src_vec;
      N_Vector rhoe_src_vec;
      N_Vector e_src_vec;
      N_Vector IR_vec;

      void *cvode_mem;
      double *dptr, *eptr, *rpar, *rparh, *abstol_ptr;
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
                        if(sundials_alloc_type%2==0)
                        {
                          if(sundials_alloc_type==0)
                            u = N_VMakeWithManagedAllocator_Cuda(neq,sunalloc,sunfree);  /* Allocate u vector */
                          else
                            u = N_VNewManaged_Cuda(neq);  /* Allocate u vector */

                          dptr=N_VGetDeviceArrayPointer_Cuda(u);

                          if(sundials_alloc_type==0)
                            e_orig = N_VMakeWithManagedAllocator_Cuda(neq,sunalloc,sunfree);  /* Allocate u vector */
                          else
                            e_orig = N_VNewManaged_Cuda(neq);  /* Allocate u vector */

                          eptr=N_VGetDeviceArrayPointer_Cuda(e_orig);
                          N_VSetCudaStream_Cuda(e_orig, &currentStream);
                          N_VSetCudaStream_Cuda(u, &currentStream);

                          if(sundials_alloc_type==0)
                            abstol_vec = N_VMakeWithManagedAllocator_Cuda(neq,sunalloc,sunfree);  
                          else
                            abstol_vec = N_VNewManaged_Cuda(neq);  /* Allocate u vector */

                          abstol_ptr = N_VGetDeviceArrayPointer_Cuda(abstol_vec);
                          N_VSetCudaStream_Cuda(abstol_vec,&currentStream);
                          amrex::Gpu::Device::streamSynchronize();
                        }
                        else
                        {
                          dptr=(double*) The_Arena()->alloc(neq*sizeof(double));
                          u = N_VMakeManaged_Cuda(neq,dptr);  /* Allocate u vector */
                          eptr= (double*) The_Arena()->alloc(neq*sizeof(double));
                          e_orig = N_VMakeManaged_Cuda(neq,eptr);  /* Allocate u vector */
                          N_VSetCudaStream_Cuda(e_orig, &currentStream);
                          N_VSetCudaStream_Cuda(u, &currentStream);

                          abstol_ptr = (double*) The_Arena()->alloc(neq*sizeof(double));
                          abstol_vec = N_VMakeManaged_Cuda(neq,abstol_ptr);
                          N_VSetCudaStream_Cuda(abstol_vec,&currentStream);
                          amrex::Gpu::streamSynchronize();
                        }
                        if(sdc_iter>=0||true)
                        {
                        T_vec = N_VMakeWithManagedAllocator_Cuda(neq,sunalloc,sunfree);
                        ne_vec = N_VMakeWithManagedAllocator_Cuda(neq,sunalloc,sunfree);
                        rho_vec = N_VMakeWithManagedAllocator_Cuda(neq,sunalloc,sunfree);
                        rho_init_vec = N_VMakeWithManagedAllocator_Cuda(neq,sunalloc,sunfree);
                        rho_src_vec = N_VMakeWithManagedAllocator_Cuda(neq,sunalloc,sunfree);
                        rhoe_src_vec = N_VMakeWithManagedAllocator_Cuda(neq,sunalloc,sunfree);
                        e_src_vec = N_VMakeWithManagedAllocator_Cuda(neq,sunalloc,sunfree);
                        IR_vec = N_VMakeWithManagedAllocator_Cuda(neq,sunalloc,sunfree);
                        }
                        amrex::Real* T_vode= N_VGetDeviceArrayPointer_Cuda(T_vec);
                        amrex::Real* ne_vode=N_VGetDeviceArrayPointer_Cuda(ne_vec);
                        amrex::Real* rho_vode=N_VGetDeviceArrayPointer_Cuda(rho_vec);
                        amrex::Real* rho_init_vode=N_VGetDeviceArrayPointer_Cuda(rho_init_vec);
                        amrex::Real* rho_src_vode=N_VGetDeviceArrayPointer_Cuda(rho_src_vec);
                        amrex::Real* rhoe_src_vode=N_VGetDeviceArrayPointer_Cuda(rhoe_src_vec);
                        amrex::Real* e_src_vode=N_VGetDeviceArrayPointer_Cuda(e_src_vec);
                        amrex::Real* IR_vode=N_VGetDeviceArrayPointer_Cuda(IR_vec);

#else
#ifdef _OPENMP
                        int nthreads=omp_get_max_threads();
                        u = N_VNew_OpenMP(neq,nthreads);  /* Allocate u vector */
                        e_orig = N_VNew_OpenMP(neq,nthreads);  /* Allocate u vector */
                        eptr=N_VGetArrayPointer_Serial(e_orig);
                        dptr=N_VGetArrayPointer_OpenMP(u);

                        abstol_vec = N_VNew_OpenMP(neq,nthreads);
                        abstol_ptr=N_VGetArrayPointer_OpenMP(abstol_vec);
                        if(sdc_iter>=0||true)
                          {
                            T_vec = N_VNew_OpenMP(neq,nthreads);
                            ne_vec = N_VNew_OpenMP(neq,nthreads);
                            rho_vec = N_VNew_OpenMP(neq,nthreads);
                            rho_init_vec = N_VNew_OpenMP(neq,nthreads);
                            rho_src_vec = N_VNew_OpenMP(neq,nthreads);
                            rhoe_src_vec = N_VNew_OpenMP(neq,nthreads);
                            e_src_vec = N_VNew_OpenMP(neq,nthreads);
                            IR_vec = N_VNew_OpenMP(neq,nthreads);
                          }
                        amrex::Real* T_vode= N_VGetArrayPointer_OpenMP(T_vec);
                        amrex::Real* ne_vode=N_VGetArrayPointer_OpenMP(ne_vec);
                        amrex::Real* rho_vode=N_VGetArrayPointer_OpenMP(rho_vec);
                        amrex::Real* rho_init_vode=N_VGetArrayPointer_OpenMP(rho_init_vec);
                        amrex::Real* rho_src_vode=N_VGetArrayPointer_OpenMP(rho_src_vec);
                        amrex::Real* rhoe_src_vode=N_VGetArrayPointer_OpenMP(rhoe_src_vec);
                        amrex::Real* e_src_vode=N_VGetArrayPointer_OpenMP(e_src_vec);
                        amrex::Real* IR_vode=N_VGetArrayPointer_OpenMP(IR_vec);
#else
                        u = N_VNew_Serial(neq);  /* Allocate u vector */
                        e_orig = N_VNew_Serial(neq);  /* Allocate u vector */
                        eptr=N_VGetArrayPointer_Serial(e_orig);
                        dptr=N_VGetArrayPointer_Serial(u);

                        abstol_vec = N_VNew_Serial(neq);
                        abstol_ptr=N_VGetArrayPointer_Serial(abstol_vec);
                        if(sdc_iter>=0||true)
                        {
                        T_vec = N_VNew_Serial(neq);
                        ne_vec = N_VNew_Serial(neq);
                        rho_vec = N_VNew_Serial(neq);
                        rho_init_vec = N_VNew_Serial(neq);
                        rho_src_vec = N_VNew_Serial(neq);
                        rhoe_src_vec = N_VNew_Serial(neq);
                        e_src_vec = N_VNew_Serial(neq);
                        IR_vec = N_VNew_Serial(neq);
                        }
                        amrex::Real* T_vode= N_VGetArrayPointer_Serial(T_vec);
                        amrex::Real* ne_vode=N_VGetArrayPointer_Serial(ne_vec);
                        amrex::Real* rho_vode=N_VGetArrayPointer_Serial(rho_vec);
                        amrex::Real* rho_init_vode=N_VGetArrayPointer_Serial(rho_init_vec);
                        amrex::Real* rho_src_vode=N_VGetArrayPointer_Serial(rho_src_vec);
                        amrex::Real* rhoe_src_vode=N_VGetArrayPointer_Serial(rhoe_src_vec);
                        amrex::Real* e_src_vode=N_VGetArrayPointer_Serial(e_src_vec);
                        amrex::Real* IR_vode=N_VGetArrayPointer_Serial(IR_vec);
#endif
#endif
      int* JH_vode_arr=NULL;
      if(inhomo_reion == 1)
          JH_vode_arr = (int*) The_Arena()->alloc(neq*sizeof(int));
      AMREX_PARALLEL_FOR_1D ( 1, i,
      {
        ode_eos_initialize_single(f_rhs_data, a, dptr, eptr, T_vode, ne_vode, rho_vode, rho_init_vode, rho_src_vode, rhoe_src_vode, e_src_vode, IR_vode, JH_vode_arr);
      });
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
                                  int idx = i+j*len.x+k*len.x*len.y-(lo.x+lo.y*len.x+lo.z*len.x*len.y);
                                  ode_eos_initialize_arrays(i, j, k, idx, f_rhs_data,
                                                            a_end, lo, len, state4, diag_eos4,
                                                            hydro_src4, reset_src4, dptr, eptr,
                                                            abstol_ptr, sdc_iter, delta_time);

#ifdef _OPENMP
                                }
                                }
                                }
#pragma omp barrier
#else
                                });
                                amrex::Gpu::Device::streamSynchronize();
#endif
                        
#ifdef CV_NEWTON
                                cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
#else
                                cvode_mem = CVodeCreate(CV_BDF);
#endif
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
                                CVodeSetUserData(cvode_mem, f_rhs_data);
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
                                //                              for (int i= 0;i < neq; ++i) {
                                  ode_eos_finalize_struct(i,j,k,idx,atomic_rates,f_rhs_data,a_end,state4,state_n4,reset_src4,diag_eos4,IR4,dptr,eptr,delta_time);
                                //PrintFinalStats(cvode_mem);
#ifdef _OPENMP
                                }
                                }
                                }
#pragma omp barrier
#else
                                });
      amrex::Gpu::Device::streamSynchronize();
#endif

    The_Arena()->free(f_rhs_data);
#ifdef AMREX_USE_CUDA
      if(sundials_alloc_type%2!=0)
      {
        The_Arena()->free(dptr);
        The_Arena()->free(eptr);
        if(use_sundials_constraint)
          The_Arena()->free(constrain);
        The_Arena()->free(abstol_ptr);
      }
#endif

                                N_VDestroy(u);          /* Free the u vector */
                                N_VDestroy(e_orig);          /* Free the e_orig vector */
                                if(use_sundials_constraint)
                                  N_VDestroy(constrain);          /* Free the constrain vector */
                                N_VDestroy(abstol_vec);          /* Free the u vector */
                                N_VDestroy(T_vec);
                                N_VDestroy(ne_vec);
                                N_VDestroy(rho_vec);
                                N_VDestroy(rho_init_vec);
                                N_VDestroy(rho_src_vec);
                                N_VDestroy(rhoe_src_vec);
                                N_VDestroy(e_src_vec);
                                N_VDestroy(IR_vec);
                                if(inhomo_reion == 1)
                                  The_Arena()->free(JH_vode_arr);
                                CVodeFree(&cvode_mem);  /* Free the integrator memory */
                              //);
                                /*                          }
                        
        }
        }*/
                                return 0;
}

#ifdef AMREX_USE_CUDA
static int f(realtype t, N_Vector u, N_Vector udot, void *user_data)
{
  Real* udot_ptr=N_VGetDeviceArrayPointer_Cuda(udot);
  Real* u_ptr=N_VGetDeviceArrayPointer_Cuda(u);
  int neq=N_VGetLength_Cuda(udot);

  auto atomic_rates = atomic_rates_glob;
  auto f_rhs_data = static_cast<RhsData*>(user_data);
  cudaStream_t currentStream = amrex::Gpu::Device::cudaStream();
  AMREX_PARALLEL_FOR_1D ( neq, idx,
  {
      f_rhs_struct(t, (u_ptr[idx]),(udot_ptr[idx]),atomic_rates,f_rhs_data,idx);
  });
  cudaStreamSynchronize(currentStream);

  return 0;
}

#else
static int f(realtype t, N_Vector u, N_Vector udot, void* user_data)
{

  Real* udot_ptr=N_VGetArrayPointer_Serial(udot);
  Real* u_ptr=N_VGetArrayPointer_Serial(u);
  int neq=N_VGetLength_Serial(udot);

  auto atomic_rates = atomic_rates_glob;
  auto f_rhs_data = static_cast<RhsData*>(user_data);
  #pragma omp parallel for
  for(int tid=0;tid<neq;tid++)
    {
                //        f_rhs_rpar(t, (u_ptr[tid]),(udot_ptr[tid]),&(rpar[4*tid]));
                f_rhs_struct(t, (u_ptr[tid]),(udot_ptr[tid]),atomic_rates,f_rhs_data,tid);
    }

  return 0;
}
#endif
