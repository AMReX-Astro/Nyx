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

using namespace amrex;

int Nyx::integrate_state_exact
  (amrex::MultiFab &S_old,
   amrex::MultiFab &D_old,
   const Real& a, const Real& delta_time)
{
    
  amrex::Gpu::setLaunchRegion(false);
  //#ifdef _OPENMP
  //#pragma omp parallel if (Gpu::notInLaunchRegion())
  //#endif
  for ( MFIter mfi(S_old, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
      double* dptr;
      //check that copy contructor vs create constructor works??
      const Box& bx = mfi.tilebox();

      const auto state = S_old.array(mfi);
      const auto diag_eos = D_old.array(mfi); 
      const Dim3 lo = amrex::lbound(bx);
      const Dim3 hi = amrex::ubound(bx);
      amrex::Print()<<lo<<hi<<bx<<std::endl;
	for (int k = lo.z; k <= hi.z; ++k) {
	  for (int j = lo.y; j <= hi.y; ++j) {
	    //	    AMREX_PRAGMA_SIMD
	      for (int i = lo.x; i <= hi.x; ++i) {
		double rparh[4];
		double e_in=state(i,j,k,Eint)/state(i,j,k,Density);
		double e_out=-10;
		rparh[0]= diag_eos(i,j,k,Temp_comp);   //rpar(1)=T_vode
		rparh[1]= diag_eos(i,j,k,Ne_comp);//    rpar(2)=ne_vode
		rparh[2]= state(i,j,k,Density); //    rpar(3)=rho_vode
		rparh[3]=1/a-1;    //    rpar(4)=z_vode
		fort_update_eos(delta_time,&e_in,&e_out,rparh);
		diag_eos(i,j,k,Temp_comp)=rparh[0];   //rpar(1)=T_vode
		diag_eos(i,j,k,Ne_comp)=rparh[1];//    rpar(2)=ne_vode
				
		state(i,j,k,Eint)  += state(i,j,k,Density) * (e_out-e_in);
		state(i,j,k,Eden)  += state(i,j,k,Density) * (e_out-e_in);

	      }
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

int Nyx::integrate_state_grownexact
  (amrex::MultiFab &S_old,
   amrex::MultiFab &D_old,
   const Real& a, const Real& delta_time)
{
    
  amrex::Gpu::setLaunchRegion(false);
  //#ifdef _OPENMP
  //#pragma omp parallel if (Gpu::notInLaunchRegion())
  //#endif
  for ( MFIter mfi(S_old, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
      double* dptr;
      //check that copy contructor vs create constructor works??
      const Box& bx = mfi.growntilebox(S_old.nGrow());

      const auto state = S_old.array(mfi);
      const auto diag_eos = D_old.array(mfi); 
      const Dim3 lo = amrex::lbound(bx);
      const Dim3 hi = amrex::ubound(bx);
      amrex::Print()<<lo<<hi<<"\n"<<bx<<std::endl;
      	for (int k = lo.z; k <= hi.z; ++k) {
	  for (int j = lo.y; j <= hi.y; ++j) {
	    //	    AMREX_PRAGMA_SIMD
	      for (int i = lo.x; i <= hi.x; ++i) {
		double rparh[4];
		double e_in=state(i,j,k,Eint)/state(i,j,k,Density);
		double e_out=-10;
		rparh[0]= diag_eos(i,j,k,Temp_comp);   //rpar(1)=T_vode
		rparh[1]= diag_eos(i,j,k,Ne_comp);//    rpar(2)=ne_vode
		rparh[2]= state(i,j,k,Density); //    rpar(3)=rho_vode
		rparh[3]=1/a-1;    //    rpar(4)=z_vode

		fort_update_eos(delta_time,&e_in,&e_out,&(rparh[0]));
		diag_eos(i,j,k,Temp_comp)=rparh[0];   //rpar(1)=T_vode
		diag_eos(i,j,k,Ne_comp)=rparh[1];//    rpar(2)=ne_vode

		state(i,j,k,Eint)  += state(i,j,k,Density) * (e_out-e_in);
		state(i,j,k,Eden)  += state(i,j,k,Density) * (e_out-e_in);

	      }
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
