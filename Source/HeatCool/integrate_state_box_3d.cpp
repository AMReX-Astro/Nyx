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

int Nyx::integrate_state_box
  (amrex::MultiFab &S_old,
   amrex::MultiFab &D_old,
   const Real& a, const Real& delta_time)
{
    // time = starting time in the simulation
  realtype reltol, abstol, t, tout, umax;
  N_Vector u;
  //  SUNLinearSolver LS;
  void *cvode_mem;
  int iout, flag;
  long int nst;

  u = NULL;
  //  LS = NULL;
  cvode_mem = NULL;

  reltol = 1e-4;  /* Set the tolerances */
  abstol = 1e-4;

  fort_ode_eos_setup(a,delta_time);
  //  amrex::Cuda::setLaunchRegion(true);
  if(S_old.nGrow()>1)
  S_old.Subtract(S_old,S_old,Eint,Eden,1,S_old.nGrow());
  else
  S_old.Subtract(S_old,S_old,Eint,Eden,1,0);
  //  amrex::Cuda::setLaunchRegion(false);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  for ( MFIter mfi(S_old, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
      double* dptr;
      t=0.0;
      //check that copy contructor vs create constructor works??
      const Box& tbx = mfi.tilebox();;
      //      if(S_old.nGrow()>1)
      //      const tbx = mfi.growntilebox(S_old.nGrow());
      
      amrex::IntVect tile_size = tbx.size();
      const int* hi = tbx.hiVect();
      const int* lo = tbx.loVect();
      int long neq2=(hi[0]-lo[0]+1-2*S_old.nGrow())*(hi[1]-lo[1]+1-2*S_old.nGrow())*(hi[2]-lo[2]+1-2*S_old.nGrow());
      //neq2=(hi[0]-lo[0]+1)*(hi[1]-lo[1]+1)*(hi[2]-lo[2]+1);
      int long neq=(tile_size[0])*(tile_size[1])*(tile_size[2]);
      
      if(neq>1)
	{
	amrex::Print()<<"Integrating a box with "<<neq<<" cells"<<std::endl;
	amrex::Print()<<"Integrating a box with "<<neq2<<"real cells"<<std::endl;
	amrex::Print()<<"Integrating a box with "<<S_old.nGrow()<<"grow cells"<<std::endl;
	amrex::Print()<<"Integrating a box with tile_size"<<tile_size<<std::endl;
	amrex::Print()<<"Integrating a box with lo tile_size"<<lo[0]<<lo[1]<<lo[2]<<std::endl;
	amrex::Print()<<"Integrating a box with hi tile_size"<<hi[0]<<hi[1]<<hi[2]<<std::endl;
	amrex::Print()<<S_old[mfi].min(Eint)<<"at index"<<S_old[mfi].minIndex(Eint)<<std::endl;
	//	neq=neq2;
	}

      /* Create a CUDA vector with initial values */
      u = N_VNew_Serial(neq);  /* Allocate u vector */
      //if(check_retval((void*)u, "N_VNew_Serial", 0)) return(1);

      /*      FSetInternalEnergy_mfab(S_old[mfi].dataPtr(),
        tbx.loVect(),
	    tbx.hiVect());  /* Initialize u vector */

      dptr=N_VGetArrayPointer_Serial(u);
      S_old[mfi].copyToMem(tbx,Nyx::Eint,1,dptr);

      /*Use N_Vector to create userdata, in order to allocate data on device*/
      N_Vector Data = N_VNew_Serial(4*neq);  // Allocate u vector 
      N_VConst(0.0,Data);
      double* rparh=N_VGetArrayPointer_Serial(Data);
      N_Vector rho_tmp = N_VNew_Serial(neq);  // Allocate u vector 
      N_VConst(0.0,rho_tmp);
      double* rho_tmp_ptr=N_VGetArrayPointer_Serial(rho_tmp);
      S_old[mfi].copyToMem(tbx,Nyx::Density,1,rho_tmp_ptr);
      N_Vector T_tmp = N_VNew_Serial(2*neq);  // Allocate u vector 
      N_VConst(0.0,T_tmp);
      double* Tne_tmp_ptr=N_VGetArrayPointer_Serial(T_tmp);
      D_old[mfi].copyToMem(tbx,Nyx::Temp_comp,2,Tne_tmp_ptr);
      
      for(int i=0;i<neq;i++)
	{
	  rparh[4*i+0]= Tne_tmp_ptr[i];   //rpar(1)=T_vode
	  rparh[4*i+1]= Tne_tmp_ptr[neq+i];//    rpar(2)=ne_vode
	  rparh[4*i+2]= rho_tmp_ptr[i]; //    rpar(3)=rho_vode
	  rparh[4*i+3]=1/a-1;    //    rpar(4)=z_vode

	}
      N_VDiv(u,rho_tmp,u);

      /* Call CVodeCreate to create the solver memory and specify the 
       * Backward Differentiation Formula and the use of a Newton iteration */
      #ifdef CV_NEWTON
      cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
      #else
      cvode_mem = CVodeCreate(CV_BDF);
      #endif
      //      if(check_retval((void *)cvode_mem, "CVodeCreate", 0)) return(1);

      /* Call CVodeInit to initialize the integrator memory and specify the
       * user's right hand side function in u'=f(t,u), the initial time T0, and
       * the initial dependent variable vector u. */
      flag = CVodeInit(cvode_mem, f, t, u);
      //      if(check_retval(&flag, "CVodeInit", 1)) return(1);

      /* Call CVodeSStolerances to specify the scalar relative tolerance
       * and scalar absolute tolerance */
      N_Vector abstol_vec = N_VNew_Serial(neq);                                                                                                                                                 
      N_VScale(abstol,u,abstol_vec);                                                                                                                                                                    
      flag = CVodeSVtolerances(cvode_mem, reltol, abstol_vec);
      //      flag = CVodeSStolerances(cvode_mem, reltol, abstol);
      //      if (check_retval(&flag, "CVodeSVtolerances", 1)) return(1);

#ifdef AMREX_USE_SUNDIALS_3x4x
      bool check_nonnegative=true;
      if(check_nonnegative)
	{
      N_Vector constrain=N_VClone(u);
      N_VConst(2,constrain);	      
      flag =CVodeSetConstraints(cvode_mem,constrain);
      //      if (check_retval(&flag, "CVodeSetConstraints", 1)) return(1);
	}
#endif

      /* Create SPGMR solver structure without preconditioning
       * and the maximum Krylov dimension maxl */
      //      LS = SUNSPGMR(u, PREC_NONE, 0);
      //      if(check_flag(&flag, "SUNSPGMR", 1)) return(1);

      /* Set CVSpils linear solver to LS */
      //      flag = CVSpilsSetLinearSolver(cvode_mem, LS);
      //      if(check_flag(&flag, "CVSpilsSetLinearSolver", 1)) return(1);

      flag = CVDiag(cvode_mem);

      CVodeSetUserData(cvode_mem, &Data);

      /* Call CVode using 1 substep */
      /*      flag = CVode(cvode_mem, tout, u, &t, CV_NORMAL);
      if(check_flag(&flag, "CVode", 1)) break;
      */
      int n_substeps=1;
      for(iout=1, tout= delta_time/n_substeps ; iout <= n_substeps; iout++, tout += delta_time/n_substeps) {
      flag = CVode(cvode_mem, tout, u, &t, CV_NORMAL);
      umax = N_VMaxNorm(u);
      flag = CVodeGetNumSteps(cvode_mem, &nst);
      check_retval(&flag, "CVodeGetNumSteps", 1);
      PrintOutput(tout, umax, nst);
      }
      
      int one_in=1;
      for(int i=0;i<neq;i++)
	{
	  Tne_tmp_ptr[i]=rparh[4*i+0];   //rpar(1)=T_vode
	  Tne_tmp_ptr[neq+i]=rparh[4*i+1];//    rpar(2)=ne_vode
	  // rho should not change  rho_tmp_ptr[i]=rparh[4*i+2]; //    rpar(3)=rho_vode
	  fort_ode_eos_finalize(&(dptr[i]), &(rparh[4*i]), one_in);
	}
      D_old[mfi].copyFromMem(tbx,Nyx::Temp_comp,2,Tne_tmp_ptr);

      N_VProd(u,rho_tmp,u);                
    
      S_old[mfi].copyFromMem(tbx,Nyx::Eint,1,dptr);
      S_old[mfi].addFromMem(tbx,Nyx::Eden,1,dptr);
      amrex::Print()<<S_old[mfi].min(Eint)<<"at index"<<S_old[mfi].minIndex(Eint)<<std::endl;
      amrex::Print()<<S_old[mfi].min(Eden)<<"at index"<<S_old[mfi].minIndex(Eden)<<std::endl;
      PrintFinalStats(cvode_mem);
      
      N_VDestroy(u);          /* Free the u vector */
      N_VDestroy(Data);          /* Free the userdata vector */
      N_VDestroy(T_tmp);
      N_VDestroy(rho_tmp);
      CVodeFree(&cvode_mem);  /* Free the integrator memory */
#ifndef NDEBUG
        if (S_old[mfi].contains_nan())
            amrex::Abort("state has NaNs after the first strang call");
#endif
    }
    return 0;
}

int Nyx::integrate_state_grownbox
  (amrex::MultiFab &S_old,
   amrex::MultiFab &D_old,
   const Real& a, const Real& delta_time)
{
    // time = starting time in the simulation
  realtype reltol, abstol, t, tout, umax;
  N_Vector u;
  //  SUNLinearSolver LS;
  void *cvode_mem;
  int iout, flag;
  long int nst;
  bool do_tiling=false;    

  u = NULL;
  //  LS = NULL;
  cvode_mem = NULL;

  reltol = 1e-4;  /* Set the tolerances */
  abstol = 1e-4;

  fort_ode_eos_setup(a,delta_time);
  if(S_old.nGrow()>1)
  S_old.Subtract(S_old,S_old,Eint,Eden,1,S_old.nGrow());
  else
  S_old.Subtract(S_old,S_old,Eint,Eden,1,0);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  for ( MFIter mfi(S_old, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
      double* dptr;
      t=0.0;
      //check that copy contructor vs create constructor works??
      //      const Box& tbx = mfi.tilebox();;
      //      if(S_old.nGrow()>1)
      const Box& tbx = mfi.growntilebox(S_old.nGrow());
      
      amrex::IntVect tile_size = tbx.size();
      const int* hi = tbx.hiVect();
      const int* lo = tbx.loVect();
      int long neq2=(hi[0]-lo[0]+1-2*S_old.nGrow())*(hi[1]-lo[1]+1-2*S_old.nGrow())*(hi[2]-lo[2]+1-2*S_old.nGrow());
      //neq2=(hi[0]-lo[0]+1)*(hi[1]-lo[1]+1)*(hi[2]-lo[2]+1);
      int long neq=(tile_size[0])*(tile_size[1])*(tile_size[2]);
      
      if(neq>1)
	{
	amrex::Print()<<"Integrating a box with "<<neq<<" cells"<<std::endl;
	amrex::Print()<<"Integrating a box with "<<neq2<<"real cells"<<std::endl;
	amrex::Print()<<"Integrating a box with "<<S_old.nGrow()<<"grow cells"<<std::endl;
	amrex::Print()<<"Integrating a box with tile_size"<<tile_size<<std::endl;
	amrex::Print()<<"Integrating a box with lo tile_size"<<lo[0]<<lo[1]<<lo[2]<<std::endl;
	amrex::Print()<<"Integrating a box with hi tile_size"<<hi[0]<<hi[1]<<hi[2]<<std::endl;
	amrex::Print()<<S_old[mfi].min(Eint)<<"at index"<<S_old[mfi].minIndex(Eint)<<std::endl;
	//	neq=neq2;
	}

      /* Create a CUDA vector with initial values */
      u = N_VNew_Serial(neq);  /* Allocate u vector */
      //      if(check_retval((void*)u, "N_VNew_Serial", 0)) return(1);

      /*      FSetInternalEnergy_mfab(S_old[mfi].dataPtr(),
        tbx.loVect(),
	    tbx.hiVect());  /* Initialize u vector */

      dptr=N_VGetArrayPointer_Serial(u);
      S_old[mfi].copyToMem(tbx,Nyx::Eint,1,dptr);

      /*Use N_Vector to create userdata, in order to allocate data on device*/
      N_Vector Data = N_VNew_Serial(4*neq);  // Allocate u vector 
      N_VConst(0.0,Data);
      double* rparh=N_VGetArrayPointer_Serial(Data);
      N_Vector rho_tmp = N_VNew_Serial(neq);  // Allocate u vector 
      N_VConst(0.0,rho_tmp);
      double* rho_tmp_ptr=N_VGetArrayPointer_Serial(rho_tmp);
      S_old[mfi].copyToMem(tbx,Nyx::Density,1,rho_tmp_ptr);
      N_Vector T_tmp = N_VNew_Serial(2*neq);  // Allocate u vector 
      N_VConst(0.0,T_tmp);
      double* Tne_tmp_ptr=N_VGetArrayPointer_Serial(T_tmp);
      D_old[mfi].copyToMem(tbx,Nyx::Temp_comp,2,Tne_tmp_ptr);
      
      for(int i=0;i<neq;i++)
	{
	  rparh[4*i+0]= Tne_tmp_ptr[i];   //rpar(1)=T_vode
	  rparh[4*i+1]= Tne_tmp_ptr[neq+i];//    rpar(2)=ne_vode
	  rparh[4*i+2]= rho_tmp_ptr[i]; //    rpar(3)=rho_vode
	  rparh[4*i+3]=1/a-1;    //    rpar(4)=z_vode

	}
      N_VDiv(u,rho_tmp,u);

      /* Call CVodeCreate to create the solver memory and specify the 
       * Backward Differentiation Formula and the use of a Newton iteration */
      #ifdef CV_NEWTON
      cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
      #else
      cvode_mem = CVodeCreate(CV_BDF);
      #endif
      //      if(check_retval((void *)cvode_mem, "CVodeCreate", 0)) return(1);

      /* Call CVodeInit to initialize the integrator memory and specify the
       * user's right hand side function in u'=f(t,u), the initial time T0, and
       * the initial dependent variable vector u. */
      flag = CVodeInit(cvode_mem, f, t, u);
      //      if(check_retval(&flag, "CVodeInit", 1)) return(1);

      /* Call CVodeSStolerances to specify the scalar relative tolerance
       * and scalar absolute tolerance */
      N_Vector abstol_vec = N_VNew_Serial(neq);                                                                                                                                                 
      N_VScale(abstol,u,abstol_vec);                                                                                                                                                                    
      flag = CVodeSVtolerances(cvode_mem, reltol, abstol_vec);
      //      flag = CVodeSStolerances(cvode_mem, reltol, abstol);
      //      if (check_retval(&flag, "CVodeSVtolerances", 1)) return(1);

#ifdef AMREX_USE_SUNDIALS_3x4x
      bool check_nonnegative=true;
      if(check_nonnegative)
	{
      N_Vector constrain=N_VClone(u);
      N_VConst(2,constrain);	      
      flag =CVodeSetConstraints(cvode_mem,constrain);
      //      if (check_retval(&flag, "CVodeSetConstraints", 1)) return(1);
	}
#endif

      /* Create SPGMR solver structure without preconditioning
       * and the maximum Krylov dimension maxl */
      //      LS = SUNSPGMR(u, PREC_NONE, 0);
      //      if(check_flag(&flag, "SUNSPGMR", 1)) return(1);

      /* Set CVSpils linear solver to LS */
      //      flag = CVSpilsSetLinearSolver(cvode_mem, LS);
      //      if(check_flag(&flag, "CVSpilsSetLinearSolver", 1)) return(1);

      flag = CVDiag(cvode_mem);

      CVodeSetUserData(cvode_mem, &Data);

      /* Call CVode using 1 substep */
      /*      flag = CVode(cvode_mem, tout, u, &t, CV_NORMAL);
      if(check_flag(&flag, "CVode", 1)) break;
      */
      int n_substeps=1;
      for(iout=1, tout= delta_time/n_substeps ; iout <= n_substeps; iout++, tout += delta_time/n_substeps) {
      flag = CVode(cvode_mem, tout, u, &t, CV_NORMAL);
      umax = N_VMaxNorm(u);
      flag = CVodeGetNumSteps(cvode_mem, &nst);
      check_retval(&flag, "CVodeGetNumSteps", 1);
      PrintOutput(tout, umax, nst);
      }
      
      int one_in=1;
      for(int i=0;i<neq;i++)
	{
	  Tne_tmp_ptr[i]=rparh[4*i+0];   //rpar(1)=T_vode
	  Tne_tmp_ptr[neq+i]=rparh[4*i+1];//    rpar(2)=ne_vode
	  // rho should not change  rho_tmp_ptr[i]=rparh[4*i+2]; //    rpar(3)=rho_vode
	  fort_ode_eos_finalize(&(dptr[i]), &(rparh[4*i]), one_in);
	}
      D_old[mfi].copyFromMem(tbx,Nyx::Temp_comp,2,Tne_tmp_ptr);

      N_VProd(u,rho_tmp,u);                
    
      S_old[mfi].copyFromMem(tbx,Nyx::Eint,1,dptr);
      S_old[mfi].addFromMem(tbx,Nyx::Eden,1,dptr);
      amrex::Print()<<S_old[mfi].min(Eint)<<"at index"<<S_old[mfi].minIndex(Eint)<<std::endl;
      amrex::Print()<<S_old[mfi].min(Eden)<<"at index"<<S_old[mfi].minIndex(Eden)<<std::endl;
      PrintFinalStats(cvode_mem);
      
      N_VDestroy(u);          /* Free the u vector */
      N_VDestroy(Data);          /* Free the userdata vector */
      N_VDestroy(T_tmp);
      N_VDestroy(rho_tmp);
      CVodeFree(&cvode_mem);  /* Free the integrator memory */
#ifndef NDEBUG
        if (S_old[mfi].contains_nan())
            amrex::Abort("state has NaNs after the first strang call");
#endif
    }
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
