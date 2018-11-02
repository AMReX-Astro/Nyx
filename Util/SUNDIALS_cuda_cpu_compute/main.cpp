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

    bool test_ic;
    int n_cell, max_grid_size;
    int cvode_meth, cvode_itmeth, write_plotfile;
    bool do_tiling;

  realtype reltol, abstol, t, tout, umax;
  N_Vector u;
  SUNLinearSolver LS;
  void *cvode_mem;
  int iout, flag;
  long int nst;

  test_ic=false;
  u = NULL;
  LS = NULL;
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

    // Create MultiFab with no ghost cells.
    MultiFab S_old(ba, dm, 7, 0);

    // Create MultiFab with no ghost cells.
    MultiFab D_old(ba, dm, 2, 0);

    
    mf.setVal(0.0);
    mf.setVal(6.226414794921875E+02);
    /*
	  rparh[4*i+0]= 3.255559960937500E+04;   //rpar(1)=T_vode
	  rparh[4*i+1]= 1.076699972152710E+00;//    rpar(2)=ne_vode
	  rparh[4*i+2]=  2.119999946752000E+12; //    rpar(3)=rho_vode
	  rparh[4*i+3]=1/(1.635780036449432E-01)-1;    //    rpar(4)=z_vode
    */
    S_old.setVal(2.119999946752000E+12,0); //    rpar(3)=rho_vode
    D_old.setVal(3.255559960937500E+04 ,0,1); //    rpar(1)=T_vode
    D_old.setVal(1.076699972152710E+00 ,1,1); //    rpar(2)=ne_vode


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
      int long neq1=(hi[0]-lo[0]+1)*(hi[1]-lo[1]+1)*(hi[2]-lo[2]+1);
      int long neq=(tile_size[0])*(tile_size[1])*(tile_size[2]);

      if(neq>1)
	amrex::Print()<<"Integrating a box with "<<neq<<" cels"<<std::endl;

      /* Create a CUDA vector with initial values */
      u = N_VNew_Cuda(neq);  /* Allocate u vector */
      if(check_retval((void*)u, "N_VNew_Cuda", 0)) return(1);

      /*      FSetInternalEnergy_mfab(mf[mfi].dataPtr(),
        tbx.loVect(),
	    tbx.hiVect());  /* Initialize u vector */

      dptr=N_VGetHostArrayPointer_Cuda(u);
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

      /* Call CVodeInit to initialize the integrator memory and specify the
       * user's right hand side function in u'=f(t,u), the initial time T0, and
       * the initial dependent variable vector u. */
      flag = CVodeInit(cvode_mem, f, t, u);
      if(check_retval(&flag, "CVodeInit", 1)) return(1);

      /* Call CVodeSStolerances to specify the scalar relative tolerance
       * and scalar absolute tolerance */
      flag = CVodeSStolerances(cvode_mem, reltol, abstol);
      if (check_retval(&flag, "CVodeSStolerances", 1)) return(1);

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
      double* rparh=N_VGetHostArrayPointer_Cuda(Data);
      N_Vector rho_tmp = N_VNew_Serial(neq);  // Allocate u vector 
      N_VConst(0.0,rho_tmp);
      double* rho_tmp_ptr=N_VGetArrayPointer_Serial(rho_tmp);
      S_old[mfi].copyToMem(tbx,0,1,rho_tmp_ptr);
      N_Vector T_tmp = N_VNew_Serial(2*neq);  // Allocate u vector 
      N_VConst(0.0,T_tmp);
      double* Tne_tmp_ptr=N_VGetArrayPointer_Serial(T_tmp);
      D_old[mfi].copyToMem(tbx,0,2,Tne_tmp_ptr);
      for(int i=0;i<neq;i++)
	{
	  rparh[4*i+0]= Tne_tmp_ptr[i];   //rpar(1)=T_vode
	  rparh[4*i+1]= Tne_tmp_ptr[neq+i];//    rpar(2)=ne_vode
	  rparh[4*i+2]= rho_tmp_ptr[i]; //    rpar(3)=rho_vode
	  rparh[4*i+3]=1/(1.635780036449432E-01)-1;    //    rpar(4)=z_vode

	}
      N_VCopyToDevice_Cuda(Data);

      CVodeSetUserData(cvode_mem, &Data);

      /*     N_Vector udot=N_VClone(u);
      f(t,u,udot,&Data);
      amrex::Print()<<"Max found: "<<N_VMaxNorm(udot)<<std::endl;
      break;*/
      /* Call CVode */
      for(iout=1, tout=8.839029760565609E-06  ; iout <= 1; iout++, tout += 8.839029760565609E-06) {
      flag = CVode(cvode_mem, tout, u, &t, CV_NORMAL);
      if(check_retval(&flag, "CVode", 1)) break;
      umax = N_VMaxNorm(u);
      flag = CVodeGetNumSteps(cvode_mem, &nst);
      check_retval(&flag, "CVodeGetNumSteps", 1);
      PrintOutput(t, umax, nst);
      }

      N_VCopyFromDevice_Cuda(u);
      //      N_VCopyFromDevice_Cuda(Data);

      for(int i=0;i<neq;i++)
	{
	  Tne_tmp_ptr[i]=rparh[4*i+0];   //rpar(1)=T_vode
	  Tne_tmp_ptr[neq+i]=rparh[4*i+1];//    rpar(2)=ne_vode
	  /*	  amrex::Print()<<"Temp"<<rparh[4*i]<<std::endl;
	  amrex::Print()<<"ne"<<rparh[4*i+1]<<std::endl;
	  amrex::Print()<<"rho"<<rparh[4*i+2]<<std::endl;
	  amrex::Print()<<"z"<<rparh[4*i+3]<<std::endl;*/
	  // rho should not change  rho_tmp_ptr[i]=rparh[4*i+2]; //    rpar(3)=rho_vode
	}
      D_old[mfi].copyFromMem(tbx,0,2,Tne_tmp_ptr);

      mf[mfi].copyFromMem(tbx,0,1,dptr);
      PrintFinalStats(cvode_mem);

      N_VDestroy(u);          /* Free the u vector */
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
      amrex::WriteSingleLevelPlotfile("PLT_DIAG_OUTPUT",
                                      D_old,
                                      {"Temp","Ne"},
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


static int f(realtype t, N_Vector u, N_Vector udot, void* user_data)
{
  N_VCopyFromDevice_Cuda(u);
  //  N_VCopyFromDevice_Cuda(*(static_cast<N_Vector*>(user_data)));
  Real* udot_ptr=N_VGetHostArrayPointer_Cuda(udot);
  Real* u_ptr=N_VGetHostArrayPointer_Cuda(u);
  int neq=N_VGetLength_Cuda(udot);
  double*  rpar=N_VGetHostArrayPointer_Cuda(*(static_cast<N_Vector*>(user_data)));
  for(int tid=0;tid<neq;tid++)
    {
      //    fprintf(stdout,"\nrpar[4*tid+0]=%g\n",rpar[4*tid]);
    RhsFnReal(t,&(u_ptr[tid]),&(udot_ptr[tid]),&(rpar[4*tid]),1);
    //    fprintf(stdout,"\nafter rpar[4*tid+0]=%g\n",rpar[4*tid]);
    }
  N_VCopyToDevice_Cuda(udot);
  //  N_VCopyToDevice_Cuda(*(static_cast<N_Vector*>(user_data)));
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
