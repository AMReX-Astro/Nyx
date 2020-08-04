#include "Nyx.H"
#include "Nyx_F.H"

using namespace amrex;

using std::string;

void
Nyx::update_state_with_sources( MultiFab& S_old, MultiFab& S_new,
                                MultiFab& ext_src_old, MultiFab& hydro_source,
                                MultiFab& grav, MultiFab& divu_cc,
                                amrex::Real dt, amrex::Real a_old, amrex::Real a_new)
{
    BL_PROFILE("Nyx::update_state_with_sources()");
    if (verbose && ParallelDescriptor::IOProcessor())
      std::cout << "Updating state with the hydro sources ... " << std::endl;
    MultiFab::RegionTag amrhydro_tag("HydroUpdate_" + std::to_string(level));
    if(verbose>1) {
        std::cout<<"hydro_src norm2(0)"<<hydro_source.norm2(0)<<std::endl;
        std::cout<<"hydro_src norm2(Eint)"<<hydro_source.norm2(Eint)<<std::endl;
        std::cout<<"hydro_src norm2(Eint)"<<hydro_source.norm2(Eden)<<std::endl;
}
    const amrex::Real a_half = 0.5 * (a_old + a_new);
    const amrex::Real a_half_inv = 1 / a_half;
    const amrex::Real a_oldsq = a_old * a_old;
    const amrex::Real a_newsq = a_new * a_new;
    const amrex::Real a_new_inv = 1.0 / a_new;
    const amrex::Real a_newsq_inv = 1.0 / a_newsq;
	/*
    int ng = 0;
    const amrex::Real a_fact[8] = {a_half_inv,a_new_inv,a_new_inv,a_new_inv,a_new_sq_inv,a_new_sq_inv,1.0,1.0};
    amrex::MultiFab::Copy(S_new, S_old, 0, 0, S_new.nComp(), ng);
    amrex::MultiFab::Saxpy(S_new, 1.0, hydro_src, 0, 0, S_new.nComp(), ng);
	*/

	////This set of dt should be used for Saxpy dt like setup
    for (amrex::MFIter mfi(S_new, amrex::TilingIfNotGPU()); mfi.isValid();
           ++mfi) {

        const amrex::Box& bx = mfi.tilebox();
        auto const& uin = S_old.array(mfi);
        auto const& uout = S_new.array(mfi);
        auto const& hydro_src = hydro_source.array(mfi);
        auto const& src = ext_src_old.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
		  for (int n = 0; n < uout.nComp(); ++n) {
			if(n==URHO)
			{
				uout(i,j,k,n) = uin(i,j,k,n) + hydro_src(i,j,k,n) 
					+ dt *  src(i,j,k,n) * a_half_inv;
			}
			else if(n>=UMX&&n<=UMZ)
			{
				uout(i,j,k,n) = a_old * uin(i,j,k,n) + hydro_src(i,j,k,n) + dt * src(i,j,k,n);
				uout(i,j,k,n) = uout(i,j,k,n) * a_new_inv;
			}
			else if(n==UEDEN)
			{
				uout(i,j,k,n) =  a_oldsq * uin(i,j,k,n) + hydro_src(i,j,k,n) + a_half * dt * src(i,j,k,n);
				uout(i,j,k,n) = uout(i,j,k,n) * a_newsq_inv;
			}
			else if(n==UEINT)
			{
				uout(i,j,k,n) =  a_oldsq*uin(i,j,k,n) + hydro_src(i,j,k,n) + a_half * dt * src(i,j,k,n) ;
				uout(i,j,k,n) = uout(i,j,k,n) * a_newsq_inv;
			}
			else
			{
				uout(i,j,k,n) = uin(i,j,k,n) + hydro_src(i,j,k,n) + dt * src(i,j,k,n) * a_half_inv;
			}
		  }
	   });
		amrex::Print()<<amrex::FArrayBox(uin)<<std::endl;
		amrex::Print()<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
		amrex::Print()<<amrex::FArrayBox(hydro_src)<<std::endl;
		amrex::Print()<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
		amrex::Print()<<amrex::FArrayBox(src)<<std::endl;
		amrex::Print()<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
		amrex::Print()<<amrex::FArrayBox(uout)<<std::endl;
		exit(0);
	}

    if(verbose>1) {
        std::cout<<"S_new norm2(0)"<<S_new.norm2(0)<<std::endl;
        std::cout<<"S_new norm2(Eint)"<<S_new.norm2(Eint)<<std::endl;
        std::cout<<"S_new norm2(Eint)"<<S_new.norm2(Eden)<<std::endl;
}

}
