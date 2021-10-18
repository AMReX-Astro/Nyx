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
#include <sundials/sundials_memory.h>
#include <sundials/sundials_config.h>

#include <nvector/nvector_serial.h>
#ifndef AMREX_USE_GPU
#ifdef _OPENMP
#include <nvector/nvector_openmp.h>
#endif
#endif
#ifdef AMREX_USE_CUDA
#include <sundials/sundials_cuda_policies.hpp>
#include <nvector/nvector_cuda.h>
#endif
#ifdef AMREX_USE_HIP
#include <sundials/sundials_hip_policies.hpp>
#include <nvector/nvector_hip.h>
#endif
#ifdef AMREX_USE_DPCPP
#include <sundials/sundials_sycl_policies.hpp>
#include <nvector/nvector_sycl.h>
#endif

#ifdef AMREX_USE_SUNDIALS_SUNMEMORY
#include <AMReX_SUNMemory.H>
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

  long int store_steps=new_max_sundials_steps;

  //#ifdef _OPENMP
  //#pragma omp parallel if (Gpu::notInLaunchRegion())
  //#endif
#ifdef _OPENMP
#ifdef AMREX_USE_GPU
  MFItInfo tiling;
  if (sundials_use_tiling)
      tiling = (MFItInfo().SetDynamic(true)).EnableTiling(sundials_tile_size);
  else
      tiling = MFItInfo().SetDynamic(true);
#pragma omp parallel
#else
  const auto tiling = (TilingIfNotGPU() && sundials_use_tiling) ? MFItInfo().EnableTiling(sundials_tile_size) : MFItInfo();
#endif
#else
  const auto tiling = (TilingIfNotGPU() && sundials_use_tiling) ? MFItInfo().EnableTiling(sundials_tile_size) : MFItInfo();
#endif
    std::string filename_inputs ="DataInputs."+std::to_string(nStep());
    std::string filename ="DataBADMAP."+std::to_string(nStep());
    std::string filename_chunk_prefix ="DataChunk."+std::to_string(nStep())+".";
    if(ParallelDescriptor::IOProcessor())
	{
	    {
            std::ofstream ofs_inputs(filename_inputs.c_str());
            ParmParse::dumpTable(ofs_inputs, true);
	    //	    ofs_inputs << "nyx.initial_a = "<<a<<std::endl;
   	    ofs_inputs << "nyx.initial_z = "<<1/a-1<<std::endl;
	    //   	    ofs_inputs << "nyx.final_a = "<<a_end<<std::endl;
       	    ofs_inputs << "nyx.final_z = "<<1/a_end-1<<std::endl;
    	    ofs_inputs << "nyx.fixed_dt = "<<delta_time<<std::endl;
       	    ofs_inputs << "nyx.hctest_filename_inputs = "<<filename_inputs<<std::endl;
	    ofs_inputs << "nyx.hctest_filename_badmap = "<<filename<<std::endl;
	    ofs_inputs << "nyx.hctest_filename_chunk = "<<filename_chunk_prefix<<std::endl;
       	    ofs_inputs << "nyx.hctest_endIndex = "<<(MFIter(S_old,tiling).length())<<std::endl;
	    }
	    /*
FArrayBox scal(Box(IntVect(AMREX_D_DECL(0,0,0)),IntVect(AMREX_D_DECL(2,0,0))),1);
	    auto scalarr=scal.array();
	    scalarr(0,0,0)=a;
	    scalarr(1,0,0)=a_end;
    	    scalarr(2,0,0)=delta_time;
	    amrex::Print()<<scal<<std::endl;
	    scal.writeOn(ofs);*/
	    {
		std::ofstream ofs(filename.c_str());
		grids.writeOn(ofs);
		dmap.writeOn(ofs);
	    }
	    /*
	    //removing const from these vars for testing purposes only
	    a=0.0;
	    a_end=0.0;
	    delta_time=0.0;*/
            std::ifstream ifs(filename.c_str());
	    //   	    scal.readFrom(ifs);
	    grids.readFrom(ifs);
	    dmap.readFrom(ifs);
	    /*
	    a=scalarr(0,0,0);
	    a_end=scalarr(1,0,0);
    	    delta_time=scalarr(2,0,0);*/
	}
    for ( MFIter mfi(S_old, tiling); mfi.isValid(); ++mfi)
    {

      //check that copy contructor vs create constructor works??
      const Box& tbx = mfi.tilebox();
      std::string filename_chunk = filename_chunk_prefix + std::to_string(mfi.index());
      std::ofstream ofs(filename_chunk.c_str());

      S_old[mfi].writeOn(ofs);
      D_old[mfi].writeOn(ofs);
      S_new[mfi].writeOn(ofs);
      hydro_src[mfi].writeOn(ofs);
      reset_src[mfi].writeOn(ofs);
      IR[mfi].writeOn(ofs);

      S_old[mfi].setVal(0);                                                                           
      D_old[mfi].setVal(0);                                                                           
      S_new[mfi].setVal(0);                                                                           
      hydro_src[mfi].setVal(0);                                                                       
      reset_src[mfi].setVal(0);                                                                       
      IR[mfi].setVal(0);

      //      string filename ="DataChunk"+std::to_string(nStep())+"."+std::to_string(mfi.index());
      std::ifstream ifs(filename_chunk.c_str());
    
      S_old[mfi].readFrom(ifs);
      D_old[mfi].readFrom(ifs);
      S_new[mfi].readFrom(ifs);
      hydro_src[mfi].readFrom(ifs);
      reset_src[mfi].readFrom(ifs);
      IR[mfi].readFrom(ifs);
      
      Array4<Real> const& state4 = S_old.array(mfi);
      Array4<Real> const& diag_eos4 = D_old.array(mfi);
      Array4<Real> const& state_n4 = S_new.array(mfi);
      Array4<Real> const& hydro_src4 = hydro_src.array(mfi);
      Array4<Real> const& reset_src4 = reset_src.array(mfi);
      Array4<Real> const& IR4 = IR.array(mfi);

      integrate_state_struct_mfin(state4,diag_eos4,state_n4,hydro_src4,reset_src4,
                                  IR4,tbx,a,a_end,delta_time,store_steps,new_max_sundials_steps,sdc_iter);
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
    BL_PROFILE_VAR("Nyx::reactions_alloc",var1);
    auto atomic_rates = atomic_rates_glob;
    auto f_rhs_data = (RhsData*)The_Arena()->alloc(sizeof(RhsData));
    Real gamma_minus_1 = gamma - 1.0;
    ode_eos_setup(f_rhs_data, gamma_minus_1, h_species);

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
      N_Vector T_vec;
      N_Vector ne_vec;
      N_Vector rho_vec;
      N_Vector rho_init_vec;
      N_Vector rho_src_vec;
      N_Vector rhoe_src_vec;
      N_Vector e_src_vec;
      N_Vector IR_vec;

      void *cvode_mem;
      realtype *dptr, *eptr, *rpar, *rparh, *abstol_ptr;
      Real *T_vode, *ne_vode,*rho_vode,*rho_init_vode,*rho_src_vode,*rhoe_src_vode,*e_src_vode,*IR_vode;
      realtype t=0.0;
                                
      u = NULL;
      e_orig = NULL;
      Data = NULL;
      abstol_vec = NULL;
      cvode_mem = NULL;

#ifdef AMREX_USE_SUNDIALS_SUNMEMORY
#ifdef AMREX_USE_DPCPP
      SUNMemoryHelper S = SUNMemoryHelper_Sycl(&amrex::Gpu::Device::streamQueue());
#endif
#endif
      long int neq = len.x*len.y*len.z;
      amrex::Gpu::streamSynchronize();
      int loop = 1;

#ifdef AMREX_USE_GPU
#ifdef AMREX_USE_CUDA
      cudaStream_t currentStream = amrex::Gpu::Device::cudaStream();
      if(sundials_alloc_type%2==0)
      {
                if(sundials_alloc_type==0)
                  u = N_VMakeWithManagedAllocator_Cuda(neq,sunalloc,sunfree);  /* Allocate u vector */
                else if(sundials_alloc_type==2)
                  u = N_VNewManaged_Cuda(neq);  /* Allocate u vector */
                else
                  u = N_VNew_Cuda(neq);  /* Allocate u vector */
                N_VSetCudaStream_Cuda(u, &currentStream);
                amrex::Gpu::Device::streamSynchronize();
      }
      else
      {
        if(sundials_alloc_type==1)
          {
                dptr=(realtype*) The_Arena()->alloc(neq*sizeof(realtype));
                u = N_VMakeManaged_Cuda(neq,dptr);  /* Allocate u vector */
                eptr= (realtype*) The_Arena()->alloc(neq*sizeof(realtype));
                e_orig = N_VMakeManaged_Cuda(neq,eptr);  /* Allocate u vector */
                N_VSetCudaStream_Cuda(e_orig, &currentStream);
                N_VSetCudaStream_Cuda(u, &currentStream);

                abstol_ptr = (realtype*) The_Arena()->alloc(neq*sizeof(realtype));
                abstol_vec = N_VMakeManaged_Cuda(neq,abstol_ptr);
                N_VSetCudaStream_Cuda(abstol_vec,&currentStream);
                T_vode=(realtype*) The_Arena()->alloc(neq*sizeof(realtype));
                T_vec = N_VMakeManaged_Cuda(neq, T_vode);
                ne_vode=(realtype*) The_Arena()->alloc(neq*sizeof(realtype));
                ne_vec = N_VMakeManaged_Cuda(neq, ne_vode);
                rho_vode=(realtype*) The_Arena()->alloc(neq*sizeof(realtype));
                rho_vec = N_VMakeManaged_Cuda(neq, rho_vode);
                rho_init_vode=(realtype*) The_Arena()->alloc(neq*sizeof(realtype));
                rho_init_vec = N_VMakeManaged_Cuda(neq, rho_init_vode);
                rho_src_vode=(realtype*) The_Arena()->alloc(neq*sizeof(realtype));
                rho_src_vec = N_VMakeManaged_Cuda(neq, rho_src_vode);
                rhoe_src_vode=(realtype*) The_Arena()->alloc(neq*sizeof(realtype));
                rhoe_src_vec = N_VMakeManaged_Cuda(neq, rhoe_src_vode);
                e_src_vode=(realtype*) The_Arena()->alloc(neq*sizeof(realtype));
                e_src_vec = N_VMakeManaged_Cuda(neq, e_src_vode);
                IR_vode=(realtype*) The_Arena()->alloc(neq*sizeof(realtype));
                IR_vec = N_VMakeManaged_Cuda(neq, IR_vode);
                N_VSetCudaStream_Cuda(T_vec,&currentStream);
                N_VSetCudaStream_Cuda(ne_vec,&currentStream);
                N_VSetCudaStream_Cuda(rho_vec,&currentStream);
                N_VSetCudaStream_Cuda(rho_init_vec,&currentStream);
                N_VSetCudaStream_Cuda(rho_src_vec,&currentStream);
                N_VSetCudaStream_Cuda(rhoe_src_vec,&currentStream);
                N_VSetCudaStream_Cuda(e_src_vec,&currentStream);
                N_VSetCudaStream_Cuda(IR_vec,&currentStream);
                amrex::Gpu::streamSynchronize();
          }
#ifdef AMREX_USE_SUNDIALS_SUNMEMORY
        else if(sundials_alloc_type==3)
          {
            u = N_VNewWithMemHelp_Cuda(neq, /*use_managed_mem=*/true, *amrex::sundials::The_SUNMemory_Helper());
            N_VSetCudaStream_Cuda(u, &currentStream);
          }
        else if(sundials_alloc_type==5)
          {
            u = N_VNewWithMemHelp_Cuda(neq, /*use_managed_mem=*/false, *amrex::sundials::The_SUNMemory_Helper());
            N_VSetCudaStream_Cuda(u, &currentStream);
          }
#endif
      }

      amrex::Gpu::streamSynchronize();
#elif defined(AMREX_USE_HIP)
      auto currentstream = amrex::Gpu::Device::gpuStream();
      //Choosing 256 here since this mimics Sundials default
      // Possibly attempt to match n_threads_and_blocks
      //     long N = ((tbx.numPts()+AMREX_GPU_NCELLS_PER_THREAD-1)/AMREX_GPU_NCELLS_PER_THREAD);
      //     SUNHipThreadDirectExecPolicy stream_exec_policy(AMREX_GPU_MAX_THREADS, currentstream);
      //     SUNHipGridStrideExecPolicy grid_exec_policy(AMREX_GPU_MAX_THREADS, std::max((N + AMREX_GPU_MAX_THREADS - 1) / AMREX_GPU_MAX_THREADS, static_cast<Long>(1)), currentstream);
      //     SUNHipBlockReduceExecPolicy reduce_exec_policy(AMREX_GPU_MAX_THREADS, std::max((N + AMREX_GPU_MAX_THREADS - 1) / AMREX_GPU_MAX_THREADS, static_cast<Long>(1)), currentstream);

      SUNHipThreadDirectExecPolicy stream_exec_policy(256, currentstream);
      SUNHipBlockReduceExecPolicy reduce_exec_policy(256, 0, currentstream);
      if(sundials_alloc_type%2==0)
      {
                if(sundials_alloc_type==0)
                  amrex::Abort("sundials_alloc_type=0 not implemented with Hip");
                else if(sundials_alloc_type==2)
                  u = N_VNewManaged_Hip(neq);  /* Allocate u vector */
                else
                  u = N_VNew_Hip(neq);  /* Allocate u vector */
                // might need a cuda analog to setting exec policy
                amrex::Gpu::Device::streamSynchronize();
      }
      else
      {
        if(sundials_alloc_type==1)
          {
                  amrex::Abort("sundials_alloc_type=1 not implemented with Hip");
          }
#ifdef AMREX_USE_SUNDIALS_SUNMEMORY
        else if(sundials_alloc_type==3)
          {
            u = N_VNewWithMemHelp_Hip(neq, /*use_managed_mem=*/true, *amrex::sundials::The_SUNMemory_Helper());
          }
        else if(sundials_alloc_type==5)
          {
            u = N_VNewWithMemHelp_Hip(neq, /*use_managed_mem=*/false, *amrex::sundials::The_SUNMemory_Helper());
          }
#endif
// might need a cuda analog to setting exec policy
      }
      N_VSetKernelExecPolicy_Hip(u, &stream_exec_policy, &reduce_exec_policy);

      amrex::Gpu::Device::streamSynchronize();
#elif defined(AMREX_USE_DPCPP)
      auto currentstream = amrex::Gpu::Device::streamQueue();
      //Choosing 256 here since this mimics Sundials default
      // Possibly attempt to match n_threads_and_blocks
      //     long N = ((tbx.numPts()+AMREX_GPU_NCELLS_PER_THREAD-1)/AMREX_GPU_NCELLS_PER_THREAD);
      //     SUNSyclThreadDirectExecPolicy stream_exec_policy(AMREX_GPU_MAX_THREADS);
      //     SUNSyclGridStrideExecPolicy grid_exec_policy(AMREX_GPU_MAX_THREADS, std::max((N + AMREX_GPU_MAX_THREADS - 1) / AMREX_GPU_MAX_THREADS, static_cast<Long>(1)));
      //     SUNSyclBlockReduceExecPolicy reduce_exec_policy(AMREX_GPU_MAX_THREADS, std::max((N + AMREX_GPU_MAX_THREADS - 1) / AMREX_GPU_MAX_THREADS, static_cast<Long>(1)));

      //Sycl version does not take a stream or queue
      SUNSyclThreadDirectExecPolicy stream_exec_policy(256);
      SUNSyclBlockReduceExecPolicy reduce_exec_policy(256, 0);
      if(sundials_alloc_type%2==0)
      {
                if(sundials_alloc_type==0)
                  amrex::Abort("sundials_alloc_type=0 not implemented with Sycl");
                else if(sundials_alloc_type==2)
                  u = N_VNewManaged_Sycl(neq, &currentstream);  /* Allocate u vector */
                else
                  u = N_VNew_Sycl(neq, &currentstream);  /* Allocate u vector */
                // might need a cuda analog to setting exec policy
                amrex::Gpu::Device::streamSynchronize();
      }
      else
      {
        if(sundials_alloc_type==1)
          {
                  amrex::Abort("sundials_alloc_type=1 not implemented with Sycl");
          }
#ifdef AMREX_USE_SUNDIALS_SUNMEMORY
        else if(sundials_alloc_type==3)
          {
              u = N_VNewWithMemHelp_Sycl(neq, /*use_managed_mem=*/true, *amrex::sundials::The_SUNMemory_Helper(), &amrex::Gpu::Device::streamQueue());
         //              u = N_VNewWithMemHelp_Sycl(neq, /*use_managed_mem=*/true, S, &amrex::Gpu::Device::streamQueue());
          }
        else if(sundials_alloc_type==5)
          {
              u = N_VNewWithMemHelp_Sycl(neq, /*use_managed_mem=*/false, *amrex::sundials::The_SUNMemory_Helper(), &amrex::Gpu::Device::streamQueue());
         //              u = N_VNewWithMemHelp_Sycl(neq, /*use_managed_mem=*/false, S, &amrex::Gpu::Device::streamQueue());
          }
#endif
// might need a cuda analog to setting exec policy
      }
      N_VSetKernelExecPolicy_Sycl(u, &stream_exec_policy, &reduce_exec_policy);

#endif
#else  /* else for ndef AMREX_USE_GPU */
#ifdef _OPENMP
              int nthreads=omp_get_max_threads();
              u = N_VNew_OpenMP(neq,nthreads);  /* Allocate u vector */
#else
              u = N_VNew_Serial(neq);  /* Allocate u vector */
#endif
#endif /* end AMREX_USE_GPU if */

              e_orig = N_VClone(u);  /* Allocate u vector */
              abstol_vec = N_VClone(u);
              if(sdc_iter>=0)
              {
                  T_vec = N_VClone(u);
                  ne_vec = N_VClone(u);
                  rho_vec = N_VClone(u);
                  rho_init_vec = N_VClone(u);
                  rho_src_vec = N_VClone(u);
                  rhoe_src_vec = N_VClone(u);
                  e_src_vec = N_VClone(u);
                  IR_vec = N_VClone(u);
              }
              else
              {
                  T_vec = N_VCloneEmpty(u);
                  ne_vec = N_VCloneEmpty(u);
                  rho_vec = N_VCloneEmpty(u);
                  rho_init_vec = N_VCloneEmpty(u);
                  rho_src_vec = N_VCloneEmpty(u);
                  rhoe_src_vec = N_VCloneEmpty(u);
                  e_src_vec = N_VCloneEmpty(u);
                  IR_vec = N_VCloneEmpty(u);
              }
              if(sundials_alloc_type!=1&&sundials_alloc_type!=7)
              {
#ifdef AMREX_USE_GPU
              eptr=N_VGetDeviceArrayPointer(e_orig);
              dptr=N_VGetDeviceArrayPointer(u);
              abstol_ptr=N_VGetDeviceArrayPointer(abstol_vec);
              T_vode= N_VGetDeviceArrayPointer(T_vec);
              ne_vode=N_VGetDeviceArrayPointer(ne_vec);
              rho_vode=N_VGetDeviceArrayPointer(rho_vec);
              rho_init_vode=N_VGetDeviceArrayPointer(rho_init_vec);
              rho_src_vode=N_VGetDeviceArrayPointer(rho_src_vec);
              rhoe_src_vode=N_VGetDeviceArrayPointer(rhoe_src_vec);
              e_src_vode=N_VGetDeviceArrayPointer(e_src_vec);
              IR_vode=N_VGetDeviceArrayPointer(IR_vec);
#else
              eptr=N_VGetArrayPointer(e_orig);
              dptr=N_VGetArrayPointer(u);
              abstol_ptr=N_VGetArrayPointer(abstol_vec);
              T_vode= N_VGetArrayPointer(T_vec);
              ne_vode=N_VGetArrayPointer(ne_vec);
              rho_vode=N_VGetArrayPointer(rho_vec);
              rho_init_vode=N_VGetArrayPointer(rho_init_vec);
              rho_src_vode=N_VGetArrayPointer(rho_src_vec);
              rhoe_src_vode=N_VGetArrayPointer(rhoe_src_vec);
              e_src_vode=N_VGetArrayPointer(e_src_vec);
              IR_vode=N_VGetArrayPointer(IR_vec);
#endif
              }

      int* JH_vode_arr=NULL;
      if(inhomo_reion == 1)
          JH_vode_arr = (int*) The_Arena()->alloc(neq*sizeof(int));
      amrex::Gpu::streamSynchronize();
      BL_PROFILE_VAR_STOP(var1);
      BL_PROFILE_VAR("Nyx::reactions_single_copy",var2);
      AMREX_PARALLEL_FOR_1D ( 1, i,
      {
        ode_eos_initialize_single(f_rhs_data, a, dptr, eptr, T_vode, ne_vode, rho_vode, rho_init_vode, rho_src_vode, rhoe_src_vode, e_src_vode, IR_vode, JH_vode_arr);
      });
      amrex::Gpu::streamSynchronize();
      BL_PROFILE_VAR_STOP(var2);
      BL_PROFILE_VAR("Nyx::reactions_cells_initialize",var3);
#ifdef AMREX_USE_GPU
      amrex::ParallelFor ( tbx, [=] AMREX_GPU_DEVICE (int i,int j,int k)
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
          ode_eos_initialize_arrays(i, j, k, idx, f_rhs_data,
                                    a_end, lo, len, state4, diag_eos4,
                                    hydro_src4, reset_src4, dptr, eptr,
                                    abstol_ptr, sdc_iter, delta_time);
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
      amrex::Gpu::streamSynchronize();
      BL_PROFILE_VAR_STOP(var3);
      BL_PROFILE_VAR("Nyx::reactions_cvsetup",var4);
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
            amrex::Gpu::streamSynchronize();
            BL_PROFILE_VAR_STOP(var4);
            BL_PROFILE_VAR("Nyx::reactions_cvode",var5);
            flag = CVode(cvode_mem, delta_time, u, &t, CV_NORMAL);
            amrex::Gpu::streamSynchronize();
            BL_PROFILE_VAR_STOP(var5);
            BL_PROFILE_VAR("Nyx::reactions_numsteps",varsteps);
            if(use_typical_steps)
              {
                long int nst=0;
                flag = CVodeGetNumSteps(cvode_mem, &nst);
                new_max_steps=std::max(nst,new_max_steps);
              }
            //                              amrex::Gpu::Device::streamSynchronize();
            //                              BL_PROFILE_VAR_STOP(cvode_timer2);
            amrex::Gpu::streamSynchronize();
            BL_PROFILE_VAR_STOP(varsteps);
            BL_PROFILE_VAR("Nyx::reactions_cells_finalize",var6);
#ifdef AMREX_USE_GPU
            AMREX_PARALLEL_FOR_3D ( tbx, i,j,k,
            {                          
#else
#ifdef _OPENMP
#pragma omp parallel for collapse(3)
      for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
#else
            AMREX_PARALLEL_FOR_3D ( tbx, i,j,k,
            {                                 
#endif
#endif
                int  idx= i + j*len.x + k*len.x*len.y - (lo.x+lo.y*len.x+lo.z*len.x*len.y);
                //                              for (int i= 0;i < neq; ++i) {
                ode_eos_finalize_struct(i,j,k,idx,atomic_rates,f_rhs_data,a_end,state4,state_n4,reset_src4,diag_eos4,IR4,dptr,eptr,delta_time);
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
    amrex::Gpu::streamSynchronize();
    BL_PROFILE_VAR_STOP(var6);
    BL_PROFILE_VAR("Nyx::reactions_free",var7);

    The_Arena()->free(f_rhs_data);

#ifdef AMREX_USE_CUDA
      if(sundials_alloc_type%2!=0&&sundials_alloc_type==1)
      {
          The_Arena()->free(dptr);
          The_Arena()->free(eptr);
          /*
          // This defaults to clone, so we don't own it
          if(use_sundials_constraint)
              The_Arena()->free(constrain);*/
          The_Arena()->free(abstol_ptr);
          The_Arena()->free(T_vode);
          The_Arena()->free(ne_vode);
          The_Arena()->free(rho_vode);
          The_Arena()->free(rho_init_vode);
          The_Arena()->free(rho_src_vode);
          The_Arena()->free(rhoe_src_vode);
          The_Arena()->free(e_src_vode);
          The_Arena()->free(IR_vode);
      }
#endif
              N_VDestroy(u);               /* Free the u vector */
              N_VDestroy(e_orig);          /* Free the e_orig vector */
              if(use_sundials_constraint)
                N_VDestroy(constrain);     /* Free the constrain vector */
              N_VDestroy(abstol_vec);      /* Free the u vector */
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
                                
    /* }
    }
    } */
    amrex::Gpu::streamSynchronize();
    BL_PROFILE_VAR_STOP(var7);
    return 0;
}

#ifdef AMREX_USE_GPU
static int f(realtype t, N_Vector u, N_Vector udot, void *user_data)
{
  amrex::Gpu::streamSynchronize();
  BL_PROFILE("Nyx::reactions_f");
  Real* udot_ptr=N_VGetDeviceArrayPointer(udot);
  Real* u_ptr=N_VGetDeviceArrayPointer(u);
  int neq=N_VGetLength(udot);

  auto atomic_rates = atomic_rates_glob;
  auto f_rhs_data = static_cast<RhsData*>(user_data);

  AMREX_PARALLEL_FOR_1D ( neq, idx,
  {
    f_rhs_struct(t, (u_ptr[idx]),(udot_ptr[idx]),atomic_rates,f_rhs_data,idx);
  });
  amrex::Gpu::streamSynchronize();

  return 0;
}

#else
static int f(realtype t, N_Vector u, N_Vector udot, void* user_data)
{
  BL_PROFILE("Nyx::reactions_f");
  Real* udot_ptr=N_VGetArrayPointer(udot);
  Real* u_ptr=N_VGetArrayPointer(u);
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
