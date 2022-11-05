#include <AMReX_Config.H>
#include <fstream>
#include <iomanip>

#include <cvode/cvode.h>               /* prototypes for CVODE fcts., consts. */
#include <cvode/cvode_diag.h>          /* access to CVDiag interface */
#include <sundials/sundials_types.h>   /* definition of type realtype */
#include <sundials/sundials_memory.h>
#include <sundials/sundials_config.h>

#include <nvector/nvector_serial.h>
#if (defined(_OPENMP) && !defined(AMREX_USE_GPU))
#include <nvector/nvector_openmp.h>
#endif
#ifdef AMREX_USE_CUDA
#include <nvector/nvector_cuda.h>
#endif
#ifdef AMREX_USE_HIP
#include <nvector/nvector_hip.h>
#endif
#ifdef AMREX_USE_DPCPP
#include <nvector/nvector_sycl.h>
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
#include <f_rhs_struct.H>

//#define MAKE_MANAGED 1
// using namespace amrex;

/* Functions Called by the Solver */
static int f(realtype t, N_Vector u, N_Vector udot, void *user_data);

static void PrintFinalStats(void *cvode_mem);
static void GetFinalStats(void *cvode_mem, N_Vector abstol_achieve, long int& nst, long int& netf, long int& nfe,
			  long int& nni, long int& ncfn, long int& nsetups, long int& nje, long int& ncfl, long int& nfeLS);

/* Private function to check function return values */
static int check_retval(void *flagvalue, const char *funcname, int opt);

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
#ifdef AMREX_USE_OMP
#ifdef AMREX_USE_GPU
  MFItInfo tiling;
  if (sundials_use_tiling)
      tiling = (MFItInfo().SetDynamic(true)).EnableTiling(sundials_tile_size);
  else
      tiling = MFItInfo().SetDynamic(true);
#else
  const auto tiling = (TilingIfNotGPU() && sundials_use_tiling) ? MFItInfo().EnableTiling(sundials_tile_size) : MFItInfo();
#endif
#else
  const auto tiling = (TilingIfNotGPU() && sundials_use_tiling) ? MFItInfo().EnableTiling(sundials_tile_size) : MFItInfo();
#endif

    amrex::ParmParse pp_nyx("nyx");
    long writeProc = ParallelDescriptor::IOProcessor() ? ParallelDescriptor::MyProc() : -1;
    int hctest_example_write = 0;
    int hctest_example_read = 0;
    int hctest_example_index = 0;
    int hctest_example_proc = 0;
    pp_nyx.query("hctest_example_write",hctest_example_write);
    if(pp_nyx.query("hctest_example_write_proc",writeProc))
        hctest_example_proc = ParallelDescriptor::MyProc()==writeProc;
    else
	hctest_example_proc = 1;

    if(hctest_example_write)
	hctest_example_index = nStep();

    pp_nyx.query("hctest_example_index",hctest_example_index);
    pp_nyx.query("hctest_example_read",hctest_example_read);
    int loc_nStep = hctest_example_index;

    //Default to these hard-coded values
    std::string filename_inputs ="hctest/inputs."+std::to_string(loc_nStep);
    std::string filename ="hctest/BADMAP."+std::to_string(loc_nStep);
    std::string filename_chunk_prefix ="hctest/Chunk."+std::to_string(loc_nStep)+".";

    int directory_overwrite = pp_nyx.query("hctest_filename_inputs",filename_inputs);
    directory_overwrite += pp_nyx.query("hctest_filename_badmap",filename);
    directory_overwrite += pp_nyx.query("hctest_filename_chunk",filename_chunk_prefix);

    if(hctest_example_read==0 && hctest_example_write!=0)
     {
      if(directory_overwrite == 0)
         amrex::UtilCreateCleanDirectory("hctest",true);
      else
        amrex::Print()<<"Using the following paths for hctest:\n"<<filename_inputs<<"\n"<<filename<<"\n"<<filename_chunk_prefix<<std::endl;
    }

    if(hctest_example_proc && hctest_example_write!=0)
    {
        sdc_writeOn(S_old,S_new, D_old, hydro_src, IR, reset_src, tiling, a, a_end, delta_time, store_steps, new_max_sundials_steps, sdc_iter, loc_nStep, filename_inputs, filename, filename_chunk_prefix);
    }
    if(hctest_example_proc && hctest_example_read!=0)
    {
        sdc_readFrom(S_old,S_new, D_old, hydro_src, IR, reset_src, tiling, a, a_end, delta_time, store_steps, new_max_sundials_steps, sdc_iter, loc_nStep, filename_inputs, filename, filename_chunk_prefix);
    }
#ifdef SAVE_REACT
	const amrex::Vector<std::string> react_in_names {"eptr-idx", "f_rhs_data-ptr-rho_init_vode-idx", "f_rhs_data-ptr-rhoe_src_vode-idx", "f_rhs_data-ptr-e_src_vode-idx", "abstol_ptr-idx", "f_rhs_data-ptr-a", "time_in"};
	const amrex::Vector<std::string> react_out_names {"dptr-idx", "f_rhs_data-ptr-rho_vode-idx", "f_rhs_data-ptr-T_vode-idx", "f_rhs_data-ptr-ne_vode-idx", "abstol_achieve_ptr-idx", "a_end", "delta_time"};
	const amrex::Vector<std::string> react_out_work_names {"nst", "netf", "nfe", "nni", "ncfn", "nsetups", "nje", "ncfl", "nfeLS"};

	MultiFab react_in(grids,dmap,react_in_names.size(),NUM_GROW);
	MultiFab react_out(grids,dmap,react_out_names.size(),NUM_GROW);
	MultiFab react_out_work(grids,dmap,react_out_work_names.size(),NUM_GROW);
#endif

#ifdef _OPENMP
#ifdef AMREX_USE_GPU
#pragma omp parallel
#endif
#endif
    for ( MFIter mfi(S_old, tiling); mfi.isValid(); ++mfi)
    {

      //check that copy contructor vs create constructor works??
      const Box& tbx = mfi.tilebox();

      Array4<Real> const& state4 = S_old.array(mfi);
      Array4<Real> const& diag_eos4 = D_old.array(mfi);
      Array4<Real> const& state_n4 = S_new.array(mfi);
      Array4<Real> const& hydro_src4 = hydro_src.array(mfi);
      Array4<Real> const& reset_src4 = reset_src.array(mfi);
      Array4<Real> const& IR4 = IR.array(mfi);
#ifdef SAVE_REACT
      Array4<Real> const& react_in_arr = react_in.array(mfi);
      Array4<Real> const& react_out_arr = react_out.array(mfi);
      Array4<Real> const& react_out_work_arr = react_out_work.array(mfi);
      integrate_state_struct_mfin(state4,diag_eos4,state_n4,hydro_src4,reset_src4,
                                  IR4,react_in_arr,react_out_arr,react_out_work_arr,
				  tbx,a,a_end,delta_time,store_steps,new_max_sundials_steps,sdc_iter);
#else
      integrate_state_struct_mfin(state4,diag_eos4,state_n4,hydro_src4,reset_src4,
                                  IR4,tbx,a,a_end,delta_time,store_steps,new_max_sundials_steps,sdc_iter);
#endif
    }
#ifdef SAVE_REACT
#ifdef NO_HYDRO
        Real cur_time = state[PhiGrav_Type].curTime();
#else
        Real cur_time = state[State_Type].curTime();
#endif
	auto plotfilename = Concatenate("plt_react_in", nStep(), 5);
        WriteSingleLevelPlotfile(plotfilename,
				 react_in, react_in_names,
				 Geom(), cur_time, nStep());
	plotfilename = Concatenate("plt_react_out", nStep(), 5);
        WriteSingleLevelPlotfile(plotfilename,
				 react_out, react_out_names,
				 Geom(), cur_time, nStep());
	plotfilename = Concatenate("plt_react_out_work", nStep(), 5);
        WriteSingleLevelPlotfile(plotfilename,
				 react_out_work, react_out_work_names,
				 Geom(), cur_time, nStep());
#endif
    return 0;
}

int Nyx::integrate_state_struct_mfin
  (amrex::Array4<Real> const& state4,
   amrex::Array4<Real> const& diag_eos4,
   amrex::Array4<Real> const& state_n4,
   amrex::Array4<Real> const& hydro_src4,
   amrex::Array4<Real> const& reset_src4,
   amrex::Array4<Real> const& IR4,
#ifdef SAVE_REACT
   amrex::Array4<Real> const& react_in_arr,
   amrex::Array4<Real> const& react_out_arr,
   amrex::Array4<Real> const& react_out_work_arr,
#endif
   const Box& tbx,
   const Real& a, const amrex::Real& a_end, const Real& delta_time,
   long int& old_max_steps, long int& new_max_steps,
   const int sdc_iter)
{
    BL_PROFILE("Nyx::integrate_state_struct_mfin");
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
      N_Vector abstol_achieve_vec;
      N_Vector T_vec;
      N_Vector ne_vec;
      N_Vector rho_vec;
      N_Vector rho_init_vec;
      N_Vector rho_src_vec;
      N_Vector rhoe_src_vec;
      N_Vector e_src_vec;
      N_Vector IR_vec;

      void *cvode_mem;
      realtype *dptr, *eptr, *rpar, *rparh, *abstol_ptr, *abstol_achieve_ptr;
      Real *T_vode, *ne_vode,*rho_vode,*rho_init_vode,*rho_src_vode,*rhoe_src_vode,*e_src_vode,*IR_vode;
      realtype t=0.0;

      u = NULL;
      e_orig = NULL;
      Data = NULL;
      abstol_vec = NULL;
      abstol_achieve_vec = NULL;
      cvode_mem = NULL;

#ifdef SUNDIALS_BUILD_WITH_PROFILING
      SUNProfiler sun_profiler = nullptr;
      SUNContext_GetProfiler(*amrex::sundials::The_Sundials_Context(),
                             &sun_profiler);
      SUNProfiler_Reset(sun_profiler);
#endif

      long int neq = len.x*len.y*len.z;
      amrex::Gpu::streamSynchronize();
      int loop = 1;

#ifdef AMREX_USE_GPU
#ifdef AMREX_USE_CUDA
      auto currentstream = amrex::Gpu::Device::gpuStream();
      //Choosing 256 here since this mimics Sundials default
      // Possibly attempt to match n_threads_and_blocks
      //     long N = ((tbx.numPts()+AMREX_GPU_NCELLS_PER_THREAD-1)/AMREX_GPU_NCELLS_PER_THREAD);
      //     SUNCudaThreadDirectExecPolicy stream_exec_policy(AMREX_GPU_MAX_THREADS, currentstream);
      //     SUNCudaGridStrideExecPolicy grid_exec_policy(AMREX_GPU_MAX_THREADS, std::max((N + AMREX_GPU_MAX_THREADS - 1) / AMREX_GPU_MAX_THREADS, static_cast<Long>(1)), currentstream);
      //     SUNCudaBlockReduceExecPolicy reduce_exec_policy(AMREX_GPU_MAX_THREADS, std::max((N + AMREX_GPU_MAX_THREADS - 1) / AMREX_GPU_MAX_THREADS, static_cast<Long>(1)), currentstream);

      SUNCudaExecPolicy* stream_exec_policy = new SUNCudaThreadDirectExecPolicy(256, currentstream);
      SUNCudaExecPolicy* reduce_exec_policy;
      if (sundials_atomic_reductions) {
        reduce_exec_policy = new SUNCudaBlockReduceAtomicExecPolicy(256, 0, currentstream);
      } else {
        reduce_exec_policy = new SUNCudaBlockReduceExecPolicy(256, 0, currentstream);
      }

      if(sundials_alloc_type%2==0)
      {
                if(sundials_alloc_type==0)
                  u = N_VNewWithMemHelp_Cuda(neq, 1, *amrex::sundials::The_SUNMemory_Helper(), *amrex::sundials::The_Sundials_Context());  /* Allocate u vector */
                else if(sundials_alloc_type==2)
                  u = N_VNewManaged_Cuda(neq, *amrex::sundials::The_Sundials_Context());  /* Allocate u vector */
                else
                  u = N_VNew_Cuda(neq, *amrex::sundials::The_Sundials_Context());  /* Allocate u vector */
                amrex::Gpu::Device::streamSynchronize();
      }
      else
      {
        if(sundials_alloc_type==1)
          {
                dptr=(realtype*) The_Arena()->alloc(neq*sizeof(realtype));
                u = N_VMakeManaged_Cuda(neq,dptr, *amrex::sundials::The_Sundials_Context());  /* Allocate u vector */
                eptr= (realtype*) The_Arena()->alloc(neq*sizeof(realtype));
                e_orig = N_VMakeManaged_Cuda(neq,eptr, *amrex::sundials::The_Sundials_Context());  /* Allocate u vector */
                N_VSetKernelExecPolicy_Cuda(e_orig, stream_exec_policy, reduce_exec_policy);
                N_VSetKernelExecPolicy_Cuda(u, stream_exec_policy, reduce_exec_policy);

                abstol_ptr = (realtype*) The_Arena()->alloc(neq*sizeof(realtype));
                abstol_vec = N_VMakeManaged_Cuda(neq,abstol_ptr, *amrex::sundials::The_Sundials_Context());
                N_VSetKernelExecPolicy_Cuda(abstol_vec, stream_exec_policy, reduce_exec_policy);
                abstol_achieve_ptr = (realtype*) The_Arena()->alloc(neq*sizeof(realtype));
                abstol_achieve_vec = N_VMakeManaged_Cuda(neq,abstol_achieve_ptr, *amrex::sundials::The_Sundials_Context());
                N_VSetKernelExecPolicy_Cuda(abstol_achieve_vec, stream_exec_policy, reduce_exec_policy);
                T_vode=(realtype*) The_Arena()->alloc(neq*sizeof(realtype));
                T_vec = N_VMakeManaged_Cuda(neq, T_vode, *amrex::sundials::The_Sundials_Context());
                ne_vode=(realtype*) The_Arena()->alloc(neq*sizeof(realtype));
                ne_vec = N_VMakeManaged_Cuda(neq, ne_vode, *amrex::sundials::The_Sundials_Context());
                rho_vode=(realtype*) The_Arena()->alloc(neq*sizeof(realtype));
                rho_vec = N_VMakeManaged_Cuda(neq, rho_vode, *amrex::sundials::The_Sundials_Context());
                rho_init_vode=(realtype*) The_Arena()->alloc(neq*sizeof(realtype));
                rho_init_vec = N_VMakeManaged_Cuda(neq, rho_init_vode, *amrex::sundials::The_Sundials_Context());
                rho_src_vode=(realtype*) The_Arena()->alloc(neq*sizeof(realtype));
                rho_src_vec = N_VMakeManaged_Cuda(neq, rho_src_vode, *amrex::sundials::The_Sundials_Context());
                rhoe_src_vode=(realtype*) The_Arena()->alloc(neq*sizeof(realtype));
                rhoe_src_vec = N_VMakeManaged_Cuda(neq, rhoe_src_vode, *amrex::sundials::The_Sundials_Context());
                e_src_vode=(realtype*) The_Arena()->alloc(neq*sizeof(realtype));
                e_src_vec = N_VMakeManaged_Cuda(neq, e_src_vode, *amrex::sundials::The_Sundials_Context());
                IR_vode=(realtype*) The_Arena()->alloc(neq*sizeof(realtype));
                IR_vec = N_VMakeManaged_Cuda(neq, IR_vode, *amrex::sundials::The_Sundials_Context());
                N_VSetKernelExecPolicy_Cuda(T_vec, stream_exec_policy, reduce_exec_policy);
                N_VSetKernelExecPolicy_Cuda(ne_vec, stream_exec_policy, reduce_exec_policy);
                N_VSetKernelExecPolicy_Cuda(rho_vec, stream_exec_policy, reduce_exec_policy);
                N_VSetKernelExecPolicy_Cuda(rho_init_vec, stream_exec_policy, reduce_exec_policy);
                N_VSetKernelExecPolicy_Cuda(rho_src_vec, stream_exec_policy, reduce_exec_policy);
                N_VSetKernelExecPolicy_Cuda(rhoe_src_vec, stream_exec_policy, reduce_exec_policy);
                N_VSetKernelExecPolicy_Cuda(e_src_vec, stream_exec_policy, reduce_exec_policy);
                N_VSetKernelExecPolicy_Cuda(IR_vec, stream_exec_policy, reduce_exec_policy);
                amrex::Gpu::streamSynchronize();
          }
        else if(sundials_alloc_type==3)
          {
            u = N_VNewWithMemHelp_Cuda(neq, /*use_managed_mem=*/true, *amrex::sundials::The_SUNMemory_Helper(), *amrex::sundials::The_Sundials_Context());
          }
        else if(sundials_alloc_type==5)
          {
            u = N_VNewWithMemHelp_Cuda(neq, /*use_managed_mem=*/false, *amrex::sundials::The_SUNMemory_Helper(), *amrex::sundials::The_Sundials_Context());
          }
      }
      N_VSetKernelExecPolicy_Cuda(u, stream_exec_policy, reduce_exec_policy);

      amrex::Gpu::Device::streamSynchronize();
#elif defined(AMREX_USE_HIP)
      auto currentstream = amrex::Gpu::Device::gpuStream();
      //Choosing 256 here since this mimics Sundials default
      // Possibly attempt to match n_threads_and_blocks
      //     long N = ((tbx.numPts()+AMREX_GPU_NCELLS_PER_THREAD-1)/AMREX_GPU_NCELLS_PER_THREAD);
      //     SUNHipThreadDirectExecPolicy stream_exec_policy(AMREX_GPU_MAX_THREADS, currentstream);
      //     SUNHipGridStrideExecPolicy grid_exec_policy(AMREX_GPU_MAX_THREADS, std::max((N + AMREX_GPU_MAX_THREADS - 1) / AMREX_GPU_MAX_THREADS, static_cast<Long>(1)), currentstream);
      //     SUNHipBlockReduceExecPolicy reduce_exec_policy(AMREX_GPU_MAX_THREADS, std::max((N + AMREX_GPU_MAX_THREADS - 1) / AMREX_GPU_MAX_THREADS, static_cast<Long>(1)), currentstream);

      SUNHipExecPolicy* stream_exec_policy = new SUNHipThreadDirectExecPolicy(256, currentstream);
      SUNHipExecPolicy* reduce_exec_policy;
      if (sundials_atomic_reductions) {
        reduce_exec_policy = new SUNHipBlockReduceAtomicExecPolicy(256, 0, currentstream);
      } else {
        reduce_exec_policy = new SUNHipBlockReduceExecPolicy(256, 0, currentstream);
      }

      if(sundials_alloc_type%2==0)
      {
                if(sundials_alloc_type==0)
                  amrex::Abort("sundials_alloc_type=0 not implemented with Hip");
                else if(sundials_alloc_type==2)
                  u = N_VNewManaged_Hip(neq, *amrex::sundials::The_Sundials_Context());  /* Allocate u vector */
                else
                  u = N_VNew_Hip(neq, *amrex::sundials::The_Sundials_Context());  /* Allocate u vector */
                // might need a cuda analog to setting exec policy
                amrex::Gpu::Device::streamSynchronize();
      }
      else
      {
        if(sundials_alloc_type==1)
          {
                  amrex::Abort("sundials_alloc_type=1 not implemented with Hip");
          }
        else if(sundials_alloc_type==3)
          {
            u = N_VNewWithMemHelp_Hip(neq, /*use_managed_mem=*/true, *amrex::sundials::The_SUNMemory_Helper(), *amrex::sundials::The_Sundials_Context());
          }
        else if(sundials_alloc_type==5)
          {
            u = N_VNewWithMemHelp_Hip(neq, /*use_managed_mem=*/false, *amrex::sundials::The_SUNMemory_Helper(), *amrex::sundials::The_Sundials_Context());
          }
// might need a cuda analog to setting exec policy
      }
      N_VSetKernelExecPolicy_Hip(u, stream_exec_policy, reduce_exec_policy);

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
      SUNSyclExecPolicy* stream_exec_policy = new SUNSyclThreadDirectExecPolicy(256);
      SUNSyclExecPolicy* reduce_exec_policy = new SUNSyclBlockReduceExecPolicy(256, 0);
      if(sundials_alloc_type%2==0)
      {
                if(sundials_alloc_type==0)
                  amrex::Abort("sundials_alloc_type=0 not implemented with Sycl");
                else if(sundials_alloc_type==2)
                  u = N_VNewManaged_Sycl(neq, &currentstream, *amrex::sundials::The_Sundials_Context());  /* Allocate u vector */
                else
                  u = N_VNew_Sycl(neq, &currentstream, *amrex::sundials::The_Sundials_Context());  /* Allocate u vector */
                // might need a cuda analog to setting exec policy
                amrex::Gpu::Device::streamSynchronize();
      }
      else
      {
        if(sundials_alloc_type==1)
          {
                  amrex::Abort("sundials_alloc_type=1 not implemented with Sycl");
          }
        else if(sundials_alloc_type==3)
          {
              u = N_VNewWithMemHelp_Sycl(neq, /*use_managed_mem=*/true, *amrex::sundials::The_SUNMemory_Helper(), &amrex::Gpu::Device::streamQueue(), *amrex::sundials::The_Sundials_Context());
         //              u = N_VNewWithMemHelp_Sycl(neq, /*use_managed_mem=*/true, S, &amrex::Gpu::Device::streamQueue(), *amrex::sundials::The_Sundials_Context());
          }
        else if(sundials_alloc_type==5)
          {
              u = N_VNewWithMemHelp_Sycl(neq, /*use_managed_mem=*/false, *amrex::sundials::The_SUNMemory_Helper(), &amrex::Gpu::Device::streamQueue(), *amrex::sundials::The_Sundials_Context());
         //              u = N_VNewWithMemHelp_Sycl(neq, /*use_managed_mem=*/false, S, &amrex::Gpu::Device::streamQueue(), *amrex::sundials::The_Sundials_Context());
          }
// might need a cuda analog to setting exec policy
      }
      N_VSetKernelExecPolicy_Sycl(u, stream_exec_policy, reduce_exec_policy);

#endif
#else  /* else for ndef AMREX_USE_GPU */
#ifdef _OPENMP
              int nthreads=omp_get_max_threads();
              u = N_VNew_OpenMP(neq,nthreads, *amrex::sundials::The_Sundials_Context());  /* Allocate u vector */
#else
              u = N_VNew_Serial(neq, *amrex::sundials::The_Sundials_Context());  /* Allocate u vector */
#endif
#endif /* end AMREX_USE_GPU if */

              e_orig = N_VClone(u);  /* Allocate u vector */
              abstol_vec = N_VClone(u);
              abstol_achieve_vec = N_VClone(u);
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
              abstol_achieve_ptr=N_VGetDeviceArrayPointer(abstol_achieve_vec);
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
              abstol_achieve_ptr=N_VGetArrayPointer(abstol_achieve_vec);
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
            if(verbose > 1)
                PrintFinalStats(cvode_mem);
#ifdef SAVE_REACT
	    long int nst, netf, nfe;
	    long int nni, ncfn;
	    long int nsetups, nje, ncfl, nfeLS;
	    GetFinalStats(cvode_mem, abstol_achieve_vec, nst, netf, nfe,
			  nni, ncfn, nsetups, nje, ncfl, nfeLS);
#endif
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
#ifdef SAVE_REACT
                ode_eos_save_react_arrays(i,j,k,idx,atomic_rates,f_rhs_data,a_end,state4,state_n4,reset_src4,diag_eos4,IR4,
					  react_in_arr, react_out_arr, react_out_work_arr,
					  dptr,eptr,abstol_ptr,abstol_achieve_ptr, sdc_iter, delta_time, 0.0, nst, netf, nfe,
					  nni, ncfn, nsetups, nje, ncfl, nfeLS);
#endif
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
          The_Arena()->free(abstol_achieve_ptr);
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
              N_VDestroy(abstol_achieve_vec);      /* Free the u vector */
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

#if defined(AMREX_USE_GPU)
              delete stream_exec_policy;
              delete reduce_exec_policy;
#endif

    /* }
    }
    } */

#ifdef SUNDIALS_BUILD_WITH_PROFILING
              SUNProfiler_Print(sun_profiler, stdout);
#endif

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

/* Get and print some final statistics */

static void PrintFinalStats(void *cvode_mem)
{
  long int nst, netf, nfe;
  long int nni, ncfn;
  long int nsetups, nje, ncfl, nfeLS;
  int retval;

  // CVODE stats
  retval = CVodeGetNumSteps(cvode_mem, &nst);
  check_retval(&retval, "CVodeGetNumSteps", 1);
  retval = CVodeGetNumErrTestFails(cvode_mem, &netf);
  check_retval(&retval, "CVodeGetNumErrTestFails", 1);
  retval = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  check_retval(&retval, "CVodeGetNumRhsEvals", 1);
  // Nonlinear solver stats
  retval = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  check_retval(&retval, "CVodeGetNumNonlinSolvIters", 1);
  retval = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  check_retval(&retval, "CVodeGetNumNonlinSolvConvFails", 1);
  // Linear solver stats
  retval = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  check_retval(&retval, "CVodeGetNumLinSolvSetups", 1);
  retval = CVDiagGetNumRhsEvals(cvode_mem, &nfeLS);
  check_retval(&retval, "CVodeGetNumLinRhsEvals", 1);

  amrex::Print() << "\nFinal Statistics..\n";
  amrex::Print() << "  Steps         = " << nst     << "\n";
  amrex::Print() << "  Err test fail = " << netf    << "\n";
  amrex::Print() << "  RHS evals     = " << nfe     << "\n";
  amrex::Print() << "  NLS iters     = " << nni     << "\n";
  amrex::Print() << "  NLS fails     = " << ncfn    << "\n";
  amrex::Print() << "  LS setups     = " << nsetups << "\n";
  amrex::Print() << "  LS RHS evals  = " << nfeLS   << "\n";

  return;
}

static void GetFinalStats(void *cvode_mem, N_Vector abstol_achieve, long int& nst, long int& netf, long int& nfe,
			  long int& nni, long int& ncfn, long int& nsetups, long int& nje, long int& ncfl, long int& nfeLS)
{

  int retval;
  // CVODE stats
  retval = CVodeGetNumSteps(cvode_mem, &nst);
  retval = CVodeGetNumErrTestFails(cvode_mem, &netf);
  retval = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  // Nonlinear solver stats
  retval = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  retval = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  // Linear solver stats
  retval = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  retval = CVDiagGetNumRhsEvals(cvode_mem, &nfeLS);
  retval = CVodeGetEstLocalErrors(cvode_mem, abstol_achieve);
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

