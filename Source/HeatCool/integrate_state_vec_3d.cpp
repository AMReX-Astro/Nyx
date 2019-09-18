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
#ifdef _OPENMP
#include <nvector/nvector_openmp.h>
#endif
#ifdef AMREX_USE_CUDA
#include <nvector/nvector_cuda.h>
#endif
#define PATCH 1
//#define MAKE_MANAGED 1
using namespace amrex;

/* Functions Called by the Solver */
static int f(realtype t, N_Vector u, N_Vector udot, void *user_data);

static void PrintFinalStats(void *cvode_mem);

static void PrintOutput(realtype t, realtype umax, long int nst);

/* Private function to check function return values */
static int check_retval(void *flagvalue, const char *funcname, int opt);

//static amrex::Arena* Managed_Arena;

void* sunalloc(size_t mem_size)
{
  return (void*) The_Managed_Arena()->alloc(mem_size);
}

void sunfree(void* ptr)
{
  The_Managed_Arena()->free(ptr);
}


int Nyx::integrate_state_vec
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
  for ( MFIter mfi(S_old, false); mfi.isValid(); ++mfi )
#else
  for ( MFIter mfi(S_old, TilingIfNotGPU()); mfi.isValid(); ++mfi )
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

int Nyx::integrate_state_vec_mfin
  (amrex::Array4<Real> const& state4,
   amrex::Array4<Real> const& diag_eos4,
   const Box& tbx,
   const Real& a, const Real& delta_time)
{

  realtype reltol, abstol;
  int flag;
    
  reltol = 1e-12;  /* Set the tolerances */
  abstol = 1e-12;

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
      amrex::Print()<<"len: "<<len<<"lo: "<<lo<<"neq: "<<neq<<std::endl;
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
#else
				N_Vector e_orig = N_VNewManaged_Cuda(neq);  /* Allocate u vector */
#endif
				double* eptr=N_VGetDeviceArrayPointer_Cuda(e_orig);
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
				double* eptr=N_VGetArrayPointer_Serial(e_orig);
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
				  rparh[4*idx+0]= diag_eos4(i,j,k,Temp_comp);   //rpar(1)=T_vode
				  rparh[4*idx+1]= diag_eos4(i,j,k,Ne_comp);//    rpar(2)=ne_vode
				  rparh[4*idx+2]= state4(i,j,k,Density); //    rpar(3)=rho_vode
				  rparh[4*idx+3]=1/a-1;    //    rpar(4)=z_vode
				  abstol_ptr[idx]= diag_eos4(i,j,k,Ne_comp)<1e-7||true ? state4(i,j,k,Eint)/state4(i,j,k,Density)*abstol : 1e4*state4(i,j,k,Eint)/state4(i,j,k,Density)*abstol ;
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

#ifdef CV_NEWTON
				cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
#else
				cvode_mem = CVodeCreate(CV_BDF);
#endif
				flag = CVodeInit(cvode_mem, f, t, u);

				N_VScale(abstol,u,abstol_vec);
				//				N_VConst(N_VMin(abstol_vec),abstol_vec);

				flag = CVodeSVtolerances(cvode_mem, reltol, abstol_vec);

				//				flag = CVodeSStolerances(cvode_mem, reltol, dptr[0]*abstol);
				flag = CVDiag(cvode_mem);

				CVodeSetMaxNumSteps(cvode_mem,2000);
				
				N_Vector constrain=N_VClone(u);
				N_VConst(2,constrain);	      
				flag =CVodeSetConstraints(cvode_mem,constrain);
				
				CVodeSetUserData(cvode_mem, &Data);
				//				CVodeSetMaxStep(cvode_mem, delta_time/10);
				//				BL_PROFILE_VAR("Nyx::strang_second_cvode",cvode_timer2);
				flag = CVode(cvode_mem, delta_time, u, &t, CV_NORMAL);
				//				amrex::Gpu::Device::streamSynchronize();
				//				BL_PROFILE_VAR_STOP(cvode_timer2);

#ifndef NDEBUG
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
      The_Managed_Arena()->free(constrain);
      The_Managed_Arena()->free(rparh);
      The_Managed_Arena()->free(abstol_ptr);
#endif
#endif
				N_VDestroy(u);          /* Free the u vector */
				N_VDestroy(e_orig);          /* Free the e_orig vector */
				N_VDestroy(constrain);          /* Free the constrain vector */
				N_VDestroy(abstol_vec);          /* Free the u vector */
				N_VDestroy(Data);          /* Free the userdata vector */
				CVodeFree(&cvode_mem);  /* Free the integrator memory */
			      //);
				/*			    }
			
	}
	}*/
				return 0;
}

int Nyx::integrate_state_grownvec
  (amrex::MultiFab &S_old,
   amrex::MultiFab &D_old,
   const Real& a, const Real& delta_time)
{

  fort_ode_eos_setup(a,delta_time);
  amrex::Gpu::LaunchSafeGuard lsg(true);

  
  const Real prev_time     = state[State_Type].prevTime();
  MultiFab S_old_tmp(S_old.boxArray(), S_old.DistributionMap(), NUM_STATE, NUM_GROW);
  MultiFab D_old_tmp(D_old.boxArray(), D_old.DistributionMap(), D_old.nComp(), NUM_GROW);

  MultiFab S_old_tmp2(S_old.boxArray(), S_old.DistributionMap(), NUM_STATE, NUM_GROW);
  MultiFab D_old_tmp2(D_old.boxArray(), D_old.DistributionMap(), D_old.nComp(), NUM_GROW);

  FillPatch(*this, S_old_tmp, NUM_GROW, prev_time, State_Type, 0, NUM_STATE);
  FillPatch(*this, D_old_tmp, NUM_GROW, prev_time, DiagEOS_Type, 0, D_old.nComp());

  FillPatch(*this, S_old_tmp2, NUM_GROW, prev_time, State_Type, 0, NUM_STATE);
  FillPatch(*this, D_old_tmp2, NUM_GROW, prev_time, DiagEOS_Type, 0, D_old.nComp());
  
  //#ifdef _OPENMP
  //#pragma omp parallel if (Gpu::notInLaunchRegion())
  //#endif
#ifdef _OPENMP
  for ( MFIter mfi(S_old, false); mfi.isValid(); ++mfi )
#else
  for ( MFIter mfi(S_old, TilingIfNotGPU()); mfi.isValid(); ++mfi )
#endif
    {

      //check that copy contructor vs create constructor works??
      const Box& tbx = mfi.growntilebox();
      Box  copy_bx = mfi.growntilebox();
      Box  copy_bbx = copy_bx.chop(0,13);
      Box  copy_bby = copy_bbx.chop(0,14);

      S_old[mfi].prefetchToDevice();
      D_old[mfi].prefetchToDevice();

      Array4<Real> const& state4 = S_old.array(mfi);
      Array4<Real> const& diag_eos4 = D_old.array(mfi);

      std::cout<<"e(24,14,19): "<<state4(24,14,19,Eint)<<std::endl;

      Array4<Real> const& state4_t = S_old_tmp.array(mfi);
      Array4<Real> const& diag_eos4_t = D_old_tmp.array(mfi);
      
      //integrate_state_vec_mfin(state4,diag_eos4,tbx,a,delta_time);
      
      integrate_state_vec_mfin(state4,diag_eos4,tbx,a,delta_time);
      /*
      integrate_state_vec_mfin(state4,diag_eos4,copy_bx,a,delta_time);
      integrate_state_vec_mfin(state4,diag_eos4,copy_bbx,a,delta_time);
      integrate_state_vec_mfin(state4,diag_eos4,copy_bby,a,delta_time);
      
            AMREX_FOR_3D(tbx,i,j,k,
		   {
		     state4_t(i,j,k,Eint)-=state4(i,j,k,Eint);
		     state4_t(i,j,k,Eden)-=state4(i,j,k,Eden);
		     diag_eos4_t(i,j,k,Temp_comp)-=diag_eos4(i,j,k,Temp_comp);
		     diag_eos4_t(i,j,k,Ne_comp)-=diag_eos4(i,j,k,Ne_comp);
		     });
      */
    }
  /*
    amrex::Real r = amrex::ReduceSum
	(D_old_tmp,D_old, NUM_GROW,
	 [=] AMREX_GPU_HOST_DEVICE (Box const& bx, FArrayBox const& fab, FArrayBox const& fab2) -> amrex::Real
	 {
	   const auto arr = fab.array();
	   const auto arr2 = fab2.array();
	   const Dim3 lo = amrex::lbound(bx);
	   const Dim3 hi = amrex::ubound(bx);
	   for (int k = lo.z; k <= hi.z; ++k) {
	     for (int j = lo.y; j <= hi.y; ++j) {
	       for (int i = lo.x; i <= hi.x; ++i) {
		 return abs((arr(i,j,k,Temp_comp)-arr2(i,j,k,Temp_comp)));
	       }
	     }
	   }
	   return 0;
	 });
      amrex::Gpu::streamSynchronize();
      amrex::Print()<<"Matching: "<<r<<std::endl;

  */

  //      amrex::Print()<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<S_old_tmp.max(Eint)<<"\t"<<S_old_tmp.max(Eden)<<"\t"<<D_old_tmp.max(Temp_comp)<<std::endl;

    return 0;
}

#ifdef AMREX_USE_CUDA
__device__ void f_rhs_test(Real t,double* u_ptr,Real* udot_ptr, Real* rpar, int neq)
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
      RhsFnReal(t,u_ptr+tid,udot_ptr+tid,rpar+4*tid,1);
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
  
  std::cout<<"t: "<<t<<std::endl;

  /*  N_VCopyFromDevice_Cuda(*(static_cast<N_Vector*>(user_data)));  
  double*  rparh=N_VGetHostArrayPointer_Cuda(*(static_cast<N_Vector*>(user_data)));
   
   fprintf(stdout,"\nt=%g \n\n",t);
   fprintf(stdout,"\nrparh[0]=%g \n\n",rparh[0]);
  fprintf(stdout,"\nrparh[1]=%g \n\n",rparh[1]);
  fprintf(stdout,"\nrparh[2]=%g \n\n",rparh[2]);
  fprintf(stdout,"\nrparh[3]=%g \n\n",rparh[3]);*/
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
  cudaStream_t currentStream = amrex::Gpu::Device::cudaStream();
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
  #pragma omp parallel for
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

