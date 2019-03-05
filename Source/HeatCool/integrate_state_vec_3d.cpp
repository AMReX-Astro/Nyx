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
  realtype reltol, abstol, t, tout, umax;
  int flag;
    
  reltol = 1e-4;  /* Set the tolerances */
  abstol = 1e-4;

  fort_ode_eos_setup(a,delta_time);
  amrex::Cuda::setLaunchRegion(false);
  //#ifdef _OPENMP
  //#pragma omp parallel if (Gpu::notInLaunchRegion())
  //#endif
  for ( MFIter mfi(S_old, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
      double* dptr;
      //check that copy contructor vs create constructor works??
      const Box& tbx = mfi.tilebox();
      //      Array4<Real> const& state = S_old.array(mfi);
      //      Array4<Real> const& diag_eos = D_old.array(mfi);
      const auto len = amrex::length(tbx);  // length of box
      const auto lo  = amrex::lbound(tbx);  // lower bound of box
      const auto state = (S_old[mfi]).view(lo);  // a view starting from lo
      const auto diag_eos = (D_old[mfi]).view(lo);  // a view starting from lo

      for       (int k = 0; k < len.z; ++k) {
	  // We know this is safe for simd on cpu.  So let's give compiler some help.
	//	AMREX_PRAGMA_SIMD
	for     (int j = 0; j < len.y; ++j) {
	                       {
				N_Vector u;
				//  SUNLinearSolver LS;
				void *cvode_mem;
				realtype t=0.0;
				
				u = NULL;
				//  LS = NULL;
				cvode_mem = NULL;
				long int neq = len.x;

				u = N_VNew_Serial(neq);  /* Allocate u vector */
				dptr=N_VGetArrayPointer_Serial(u);
				//if(check_retval((void*)u, "N_VNew_Serial", 0)) return(1);

				for (int i = 0; i < neq; ++i) {
				  dptr[i]=state(i,j,k,Eint)/state(i,j,k,Density);
				  state(i,j,k,Eint)  -= state(i,j,k,Density) * (dptr[i]);
				  state(i,j,k,Eden)  -= state(i,j,k,Density) * (dptr[i]);
				}
				
				N_Vector Data = N_VNew_Serial(4*neq);  // Allocate u vector 
				N_VConst(0.0,Data);
				double* rparh=N_VGetArrayPointer_Serial(Data);
				for (int i= 0;i < neq; ++i) {
				  rparh[4*i+0]= diag_eos(i,j,k,Temp_comp);   //rpar(1)=T_vode
				  rparh[4*i+1]= diag_eos(i,j,k,Ne_comp);//    rpar(2)=ne_vode
				  rparh[4*i+2]= state(i,j,k,Density); //    rpar(3)=rho_vode
				  rparh[4*i+3]=1/a-1;    //    rpar(4)=z_vode
				}
#ifdef CV_NEWTON
				cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
#else
				cvode_mem = CVodeCreate(CV_BDF);
#endif
				flag = CVodeInit(cvode_mem, f, t, u);
				N_Vector abstol_vec = N_VNew_Serial(neq);
				N_VScale(abstol,u,abstol_vec);
				flag = CVodeSVtolerances(cvode_mem, reltol, abstol_vec);
      
				//				flag = CVodeSStolerances(cvode_mem, reltol, dptr[0]*abstol);
				flag = CVDiag(cvode_mem);

				CVodeSetUserData(cvode_mem, &Data);
				flag = CVode(cvode_mem, delta_time, u, &t, CV_NORMAL);

				//				diag_eos(i,j,k,Temp_comp)=rparh[0];   //rpar(1)=T_vode
				//	diag_eos(i,j,k,Ne_comp)=rparh[1];//    rpar(2)=ne_vode
				// rho should not change  rho_tmp_ptr[i]=rparh[4*i+2]; //    rpar(3)=rho_vode

				for (int i= 0;i < neq; ++i) {
				  fort_ode_eos_finalize(&(dptr[i]), &(rparh[4*i]), 1);
				  diag_eos(i,j,k,Temp_comp)=rparh[4*i+0];   //rpar(1)=T_vode
				  diag_eos(i,j,k,Ne_comp)=rparh[4*i+1];//    rpar(2)=ne_vode
				
				  state(i,j,k,Eint)  += state(i,j,k,Density) * (dptr[i]);
				  state(i,j,k,Eden)  += state(i,j,k,Density) * (dptr[i]);
				}
				//PrintFinalStats(cvode_mem);
				
				N_VDestroy(u);          /* Free the u vector */
				N_VDestroy(abstol_vec);          /* Free the u vector */
				N_VDestroy(Data);          /* Free the userdata vector */
				CVodeFree(&cvode_mem);  /* Free the integrator memory */
			      }//);
			
	}
      }

    }
  #ifdef NDEBUG
  #ifndef NDEBUG
        if (S_old.contains_nan())
            amrex::Abort("state has NaNs after the second strang call");
  #endif
#endif

    return 0;
}

int Nyx::integrate_state_grownvec
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

  fort_ode_eos_setup(a,delta_time);
  amrex::Cuda::setLaunchRegion(false);
  //#ifdef _OPENMP
  //#pragma omp parallel if (Gpu::notInLaunchRegion())
  //#endif
  for ( MFIter mfi(S_old, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
      double* dptr;
      //check that copy contructor vs create constructor works??
      const Box& tbx = mfi.growntilebox();
            //      Array4<Real> const& state = S_old.array(mfi);
      //      Array4<Real> const& diag_eos = D_old.array(mfi);
      const auto len = amrex::length(tbx);  // length of box
      const auto lo  = amrex::lbound(tbx);  // lower bound of box
      const auto state = (S_old[mfi]).view(lo);  // a view starting from lo
      const auto diag_eos = (D_old[mfi]).view(lo);  // a view starting from lo

      for       (int k = 0; k < len.z; ++k) {
	  // We know this is safe for simd on cpu.  So let's give compiler some help.
	//	AMREX_PRAGMA_SIMD
	for     (int j = 0; j < len.y; ++j) {
	                       {
				N_Vector u;
				//  SUNLinearSolver LS;
				void *cvode_mem;
				realtype t=0.0;
				
				u = NULL;
				//  LS = NULL;
				cvode_mem = NULL;
				long int neq = len.x;

				u = N_VNew_Serial(neq);  /* Allocate u vector */
				dptr=N_VGetArrayPointer_Serial(u);
				//if(check_retval((void*)u, "N_VNew_Serial", 0)) return(1);

				for (int i = 0; i < neq; ++i) {
				  dptr[i]=state(i,j,k,Eint)/state(i,j,k,Density);
				  state(i,j,k,Eint)  -= state(i,j,k,Density) * (dptr[i]);
				  state(i,j,k,Eden)  -= state(i,j,k,Density) * (dptr[i]);
				}
				
				N_Vector Data = N_VNew_Serial(4*neq);  // Allocate u vector 
				N_VConst(0.0,Data);
				double* rparh=N_VGetArrayPointer_Serial(Data);
				for (int i= 0;i < neq; ++i) {
				  rparh[4*i+0]= diag_eos(i,j,k,Temp_comp);   //rpar(1)=T_vode
				  rparh[4*i+1]= diag_eos(i,j,k,Ne_comp);//    rpar(2)=ne_vode
				  rparh[4*i+2]= state(i,j,k,Density); //    rpar(3)=rho_vode
				  rparh[4*i+3]=1/a-1;    //    rpar(4)=z_vode
				}
#ifdef CV_NEWTON
				cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
#else
				cvode_mem = CVodeCreate(CV_BDF);
#endif
				flag = CVodeInit(cvode_mem, f, t, u);
				N_Vector abstol_vec = N_VNew_Serial(neq);
				N_VScale(abstol,u,abstol_vec);
				flag = CVodeSVtolerances(cvode_mem, reltol, abstol_vec);
      
				//				flag = CVodeSStolerances(cvode_mem, reltol, dptr[0]*abstol);
				flag = CVDiag(cvode_mem);

				CVodeSetUserData(cvode_mem, &Data);
				flag = CVode(cvode_mem, delta_time, u, &t, CV_NORMAL);

				//				diag_eos(i,j,k,Temp_comp)=rparh[0];   //rpar(1)=T_vode
				//	diag_eos(i,j,k,Ne_comp)=rparh[1];//    rpar(2)=ne_vode
				// rho should not change  rho_tmp_ptr[i]=rparh[4*i+2]; //    rpar(3)=rho_vode

				for (int i= 0;i < neq; ++i) {
				  fort_ode_eos_finalize(&(dptr[i]), &(rparh[4*i]), 1);
				  diag_eos(i,j,k,Temp_comp)=rparh[4*i+0];   //rpar(1)=T_vode
				  diag_eos(i,j,k,Ne_comp)=rparh[4*i+1];//    rpar(2)=ne_vode
				
				  state(i,j,k,Eint)  += state(i,j,k,Density) * (dptr[i]);
				  state(i,j,k,Eden)  += state(i,j,k,Density) * (dptr[i]);
				}
				//PrintFinalStats(cvode_mem);
				
				N_VDestroy(u);          /* Free the u vector */
				N_VDestroy(abstol_vec);          /* Free the u vector */
				N_VDestroy(Data);          /* Free the userdata vector */
				CVodeFree(&cvode_mem);  /* Free the integrator memory */
			      }//);
			
	}
      }
    }
  #ifndef NDEBUG
        if (S_old.contains_nan())
            amrex::Abort("state has NaNs after the first strang call");
  #endif

    return 0;
}

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
