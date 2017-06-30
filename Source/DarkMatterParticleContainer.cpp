#include "DarkMatterParticleContainer.H"
#include "dm_F.H"

using namespace amrex;

void
DarkMatterParticleContainer::moveKickDrift (amrex::MultiFab&       acceleration,
		                            int                    lev,
                    			    amrex::Real            dt,
		                	    amrex::Real            a_old,
					    amrex::Real            a_half,
					    int                    where_width)
{
    BL_PROFILE("DarkMatterParticleContainer::moveKickDrift()");

    //If there are no particles at this level
    if (lev >= this->GetParticles().size())
        return;

    const Real* dx = Geom(lev).CellSize();

    amrex::MultiFab* ac_ptr;
    if (this->OnSameGrids(lev, acceleration))
    {
        ac_ptr = &acceleration;
    }
    else
    {
        ac_ptr = new amrex::MultiFab(this->m_gdb->ParticleBoxArray(lev),
			             this->m_gdb->ParticleDistributionMap(lev),
				     acceleration.nComp(),acceleration.nGrow());
        for (amrex::MFIter mfi(*ac_ptr); mfi.isValid(); ++mfi)
            ac_ptr->setVal(0.);
        ac_ptr->copy(acceleration,0,0,acceleration.nComp());
        ac_ptr->FillBoundary();
    }

    const Real* plo = Geom(lev).ProbLo();

    int do_move = 1;

    for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {

        AoS& particles = pti.GetArrayOfStructs();
        int Np = particles.size();

        if (Np > 0)
        {
           const Box& ac_box = (*ac_ptr)[pti].box();

           update_dm_particles(&Np, particles.data(),
                               (*ac_ptr)[pti].dataPtr(),
                               ac_box.loVect(), ac_box.hiVect(),
                               plo,dx,dt,a_old,a_half,&do_move);
        }
    }

    if (ac_ptr != &acceleration) delete ac_ptr;
    
    ParticleLevel&    pmap          = this->GetParticles(lev);
    if (lev > 0 && sub_cycle)
    {
        amrex::ParticleLocData pld; 
        for (auto& kv : pmap) {
            AoS&  pbox       = kv.second.GetArrayOfStructs();
            const int   n    = pbox.size();

#ifdef _OPENMP
#pragma omp parallel for private(pld)
#endif
            for (int i = 0; i < n; i++)
            {
                ParticleType& p = pbox[i];
                if (p.id() <= 0) continue;

                // Move the particle to the proper ghost cell. 
                //      and remove any *ghost* particles that have gone too far
                // Note that this should only negate ghost particles, not real particles.
                if (!this->Where(p, pld, lev, lev, where_width))
                {
                    // Assert that the particle being removed is a ghost particle;
                    // the ghost particle is no longer in relevant ghost cells for this grid.
                    if (p.id() == amrex::GhostParticleID)
                    {
                        p.id() = -1;
                    }
                    else
                    {
                        std::cout << "Oops -- removing particle " << p.id() << std::endl;
                        amrex::Error("Trying to get rid of a non-ghost particle in moveKickDrift");
                    }
                }
            }
        }
    }
}

void
DarkMatterParticleContainer::moveKick (MultiFab&       acceleration,
                                       int             lev,
                                       Real            dt,
                                       Real            a_new,
                                       Real            a_half) 
{
    BL_PROFILE("DarkMatterParticleContainer::moveKick()");

    const Real* dx = Geom(lev).CellSize();

    MultiFab* ac_ptr;
    if (OnSameGrids(lev,acceleration))
    {
        ac_ptr = &acceleration;
    }
    else 
    {
        ac_ptr = new MultiFab(ParticleBoxArray(lev),
				  ParticleDistributionMap(lev),
				  acceleration.nComp(),acceleration.nGrow());
        for (MFIter mfi(*ac_ptr); mfi.isValid(); ++mfi)
            ac_ptr->setVal(0.);
        ac_ptr->copy(acceleration,0,0,acceleration.nComp());
        ac_ptr->FillBoundary();
    }

    const Real* plo = Geom(lev).ProbLo();

    int do_move = 0;

    for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {

        AoS& particles = pti.GetArrayOfStructs();
        int Np = particles.size();

        if (Np > 0)
        {
           const Box& ac_box = (*ac_ptr)[pti].box();

           update_dm_particles(&Np, particles.data(),
                               (*ac_ptr)[pti].dataPtr(),
                               ac_box.loVect(), ac_box.hiVect(),
                               plo,dx,dt,a_half,a_new,&do_move);
        }
    }
    
    if (ac_ptr != &acceleration) delete ac_ptr;
}

void
DarkMatterParticleContainer::InitCosmo1ppcMultiLevel(
                        MultiFab& mf, const Real disp_fac[], const Real vel_fac[], 
                        const Real particleMass, int disp_idx, int vel_idx, 
                        BoxArray &baWhereNot, int lev, int nlevs)
{
    BL_PROFILE("DarkMatterParticleContainer::InitCosmo1ppcMultiLevel()");
    const int       MyProc   = ParallelDescriptor::MyProc();
    const Geometry& geom     = m_gdb->Geom(lev);
    const Real*     dx       = geom.CellSize();

    static Array<int> calls;

    calls.resize(nlevs);

    calls[lev]++;

    if (calls[lev] > 1) return;

    Array<ParticleLevel>& particles = this->GetParticles();

    particles.reserve(15);  // So we don't ever have to do any copying on a resize.

    particles.resize(nlevs);

    ParticleType p;
    Real         disp[BL_SPACEDIM];
    Real         vel[BL_SPACEDIM];
    
    Real 	mean_disp[BL_SPACEDIM]={D_DECL(0,0,0)};


    //
    // The mf should be initialized according to the ics...
    //
    int outside_counter=0;
    long outcount[3]={0,0,0};
    long outcountminus[3]={0,0,0};
    long totalcount=0;
    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        FArrayBox&  myFab  = mf[mfi];
	const Box&  vbx    = mfi.validbox();
        const int  *fab_lo = vbx.loVect();
        const int  *fab_hi = vbx.hiVect();
        ParticleLocData pld;
        for (int kx = fab_lo[2]; kx <= fab_hi[2]; kx++)
        {
            for (int jx = fab_lo[1]; jx <= fab_hi[1]; jx++)
            {
                for (int ix = fab_lo[0]; ix <= fab_hi[0]; ix++)
                {
            	    IntVect indices(D_DECL(ix, jx, kx));
		    totalcount++;
		    if (baWhereNot.contains(indices)) 
		    {
                       continue;
		    }

	            for (int n = 0; n < BL_SPACEDIM; n++)
	            {
                        disp[n] = myFab(indices,disp_idx+n);
                        //
			// Start with homogeneous distribution (for 1 p per cell in the center of the cell),
			//
	                p.pos(n) = geom.ProbLo(n) + 
                            (indices[n]+Real(0.5))*dx[n];
			if(disp[n]*disp_fac[n]>dx[n]/2.0)
			  outcount[n]++;
			if(disp[n]*disp_fac[n]<-dx[n]/2.0)
			  outcountminus[n]++;
			mean_disp[n]+=fabs(disp[n]);
			//
                        // then add the displacement (input values weighted by domain length).
                        //
	                p.pos(n) += disp[n] * disp_fac[n];

                        //
			// Set the velocities.
                        //
                        vel[n] = myFab(indices,vel_idx+n);
	                p.rdata(n+1) = vel[n] * vel_fac[n];
	            }
                    //
		    // Set the mass of the particle from the input value.
                    //
	            p.rdata(0)  = particleMass;
	            p.id()      = ParticleType::NextID();
	            p.cpu()     = MyProc;
	
	            if (!this->Where(p, pld))
                    {
      		        this->PeriodicShift(p);

                        if (!this->Where(p, pld))
                            amrex::Abort("DarkMatterParticleContainer::InitCosmo1ppcMultiLevel():invalid particle");
                    }

		    BL_ASSERT(pld.m_lev >= 0 && pld.m_lev <= m_gdb->finestLevel());
		    //handle particles that ran out of this level into a finer one. 
		    if (baWhereNot.contains(pld.m_cell))
		    {
		      outside_counter++;
		      ParticleType newp[8];
                      ParticleLocData new_pld;
		      for (int i=0;i<8;i++)
		      {
                          newp[i].rdata(0)   = particleMass/8.0;
                          newp[i].id()       = ParticleType::NextID();
                          newp[i].cpu()      = MyProc;
                          for (int dim=0;dim<BL_SPACEDIM;dim++)
                          {
                              newp[i].pos(dim)=p.pos(dim)+(2*((i/(1 << dim)) % 2)-1)*dx[dim]/4.0;
                              newp[i].rdata(dim+1)=p.rdata(dim+1);
                          }
                          
                          if (!this->Where(newp[i], new_pld))
                          {
                              this->PeriodicShift(newp[i]);
                              
                              if (!this->Where(newp[i], new_pld))
                                  amrex::Abort("DarkMatterParticleContainer::InitCosmo1ppcMultiLevel():invalid particle");
                          }
                          particles[new_pld.m_lev][std::make_pair(new_pld.m_grid, 
                                                                  new_pld.m_tile)].push_back(newp[i]);
		      }
		      
		    }
	            
	            //
	            // Add it to the appropriate PBox at the appropriate level.
	            //
		    else
                        particles[pld.m_lev][std::make_pair(pld.m_grid, pld.m_tile)].push_back(p);
                }
            }
        }
    }
    Redistribute();
}

void
DarkMatterParticleContainer::InitCosmo1ppc(MultiFab& mf, const Real vel_fac[], const Real particleMass)
{
    BL_PROFILE("DarkMatterParticleContainer::InitCosmo1ppc()");
    const int       MyProc   = ParallelDescriptor::MyProc();
    const Geometry& geom     = m_gdb->Geom(0);
    const Real*     dx       = geom.CellSize();

    Array<ParticleLevel>& particles = this->GetParticles();

    particles.reserve(15);  // So we don't ever have to do any copying on a resize.

    particles.resize(m_gdb->finestLevel()+1);

    for (int lev = 0; lev < particles.size(); lev++)
    {
        BL_ASSERT(particles[lev].empty());
    }

    ParticleType      p;
    ParticleLocData   pld;
    Real              disp[BL_SPACEDIM];
    const Real        len[BL_SPACEDIM] = { D_DECL(geom.ProbLength(0),
                                                  geom.ProbLength(1),
                                                  geom.ProbLength(2)) };
    //
    // The grid should be initialized according to the ics...
    //
    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        FArrayBox&  myFab  = mf[mfi];
	const Box&  vbx    = mfi.validbox();
        const int  *fab_lo = vbx.loVect();
        const int  *fab_hi = vbx.hiVect();

        for (int kx = fab_lo[2]; kx <= fab_hi[2]; kx++)
        {
            for (int jx = fab_lo[1]; jx <= fab_hi[1]; jx++)
            {
                for (int ix = fab_lo[0]; ix <= fab_hi[0]; ix++)
                {
            	    IntVect indices(D_DECL(ix, jx, kx));

	            for (int n = 0; n < BL_SPACEDIM; n++)
	            {
                        disp[n] = myFab(indices,n);
                        //
			// Start with homogeneous distribution (for 1 p per cell in the center of the cell),
                        // then add the displacement (input values weighted by domain length).
                        //
	                p.pos(n) = geom.ProbLo(n) + 
                            (indices[n]+Real(0.5))*dx[n] +
                            disp[n] * len[n];
                        //
			// Set the velocities.
                        //
	                p.rdata(n+1) = disp[n] * vel_fac[n];
	            }
                    //
		    // Set the mass of the particle from the input value.
                    //
	            p.rdata(0)  = particleMass;
	            p.id()      = ParticleType::NextID();
	            p.cpu()     = MyProc;
	
	            if (!this->Where(p, pld))
                    {
      		        this->PeriodicShift(p);
                        
                        if (!this->Where(p, pld))
                            amrex::Abort("DarkMatterParticleContainer::InitCosmo1ppc(): invalid particle");
		    }

	            BL_ASSERT(pld.m_lev >= 0 && pld.m_lev <= this->finestLevel());
	            //
	            // Add it to the appropriate PBox at the appropriate level.
	            //
	            particles[pld.m_lev][std::make_pair(pld.m_grid, pld.m_tile)].push_back(p);
                }
            }
        }
    }
}

void
DarkMatterParticleContainer::InitCosmo(
            MultiFab& mf, const Real vel_fac[], const Array<int> n_part, const Real particleMass)
{
    Real shift[] = {0,0,0};
    InitCosmo(mf, vel_fac, n_part, particleMass, shift);
}

void
DarkMatterParticleContainer::InitCosmo(
            MultiFab& mf, const Real vel_fac[], const Array<int> n_part, const Real particleMass, const Real shift[])
{
    BL_PROFILE("DarkMatterParticleContainer::InitCosmo()");
    const int       MyProc   = ParallelDescriptor::MyProc();
    const int       IOProc   = ParallelDescriptor::IOProcessorNumber();
    const Real      strttime = ParallelDescriptor::second();
    const Geometry& geom     = m_gdb->Geom(0);

    Array<ParticleLevel>& particles = this->GetParticles();

    particles.reserve(15);  // So we don't ever have to do any copying on a resize.

    particles.resize(m_gdb->finestLevel()+1);

    for (int lev = 0; lev < particles.size(); lev++)
    {
        BL_ASSERT(particles[lev].empty());
    }

    const Real len[BL_SPACEDIM] = { D_DECL(geom.ProbLength(0),
                                           geom.ProbLength(1),
                                           geom.ProbLength(2)) };
    //
    // Print the grids as a sanity check.
    //
    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
	const Box&  vbx    = mfi.validbox();
        const int  *fab_lo = vbx.loVect();
        const int  *fab_hi = vbx.hiVect();
        if (vbx.isEmpty())
        {
           std::cout << "...bad grid lo " << fab_lo[0] << " " << fab_lo[1] << " " << fab_lo[2] << '\n';
           std::cout << "...bad grid hi " << fab_hi[0] << " " << fab_hi[1] << " " << fab_hi[2] << '\n';
           amrex::Error("Empty box in InitCosmo ");
        }
        if (!geom.Domain().contains(vbx))
        {
           std::cout << "...bad grid lo " << fab_lo[0] << " " << fab_lo[1] << " " << fab_lo[2] << '\n';
           std::cout << "...bad grid hi " << fab_hi[0] << " " << fab_hi[1] << " " << fab_hi[2] << '\n';
           amrex::Error("Box in InitCosmo not contained in domain");
        }
    }
    //
    // We will need one ghost cell, so check wether we have one.
    //
    if (mf.nGrow() < 1)
        amrex::Abort("DarkMatterParticleContainer::InitCosmo: mf needs at least one correctly filled ghost zone!");

    if ( !(n_part[0] == n_part[1] && n_part[1] == n_part[2]) )
    {
	    std::cout << '\n' << '\n';
	    std::cout << "Your particle lattice will have different spacings in the spatial directions!" << '\n';
	    std::cout << "You might want to change the particle number or the algorithm... ;)" << '\n';
	    std::cout << '\n' << '\n';
    }
    //
    // Place the particles evenly spaced in the problem domain.
    // Not perfectly fast - but easy
    //
    Real         pos[BL_SPACEDIM];
    ParticleType p;

    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        const Box&  box     = mfi.validbox();
        RealBox     gridloc = RealBox(box, geom.CellSize(), geom.ProbLo());
	const Real* xlo     = gridloc.lo();
	const Real* xhi     = gridloc.hi();

        ParticleLocData pld;
        for (int k = 0; k < n_part[2]; k++)
        {
	    for (int j = 0; j < n_part[1]; j++)
            {
	        for (int i = 0; i < n_part[0]; i++)
                {
		    bool    isInValidBox = true;
            	    IntVect indices(D_DECL(i, j, k));

		    for (int n = 0; n < BL_SPACEDIM; n++)
                    {
  		        pos[n] = geom.ProbLo(n)
		               + (indices[n] + Real(0.5))*len[n]/n_part[n]
			       + shift[n];
                        //
			// Make sure particle is not on a boundary...
                        //
			pos[n] += 1e-14 * (geom.ProbHi(n) - geom.ProbLo(n));

			isInValidBox = isInValidBox 
				     && (pos[n] > xlo[n]) 
				     && (pos[n] < xhi[n]);
		    }

		    if (isInValidBox)
                    {
                        D_TERM(p.pos(0) = pos[0];,
                               p.pos(1) = pos[1];,
                               p.pos(2) = pos[2];);
                        //
		        // Set the mass of the particle.
                        //
	                p.rdata(0)  = particleMass;
	                p.id()      = ParticleType::NextID();
	                p.cpu()     = MyProc;

	                if (!this->Where(p, pld))
                        {
      		            this->PeriodicShift(p);

                            if (!this->Where(p, pld))
                                amrex::Abort("DarkMatterParticleContainer::InitCosmo(): invalid particle");
		        }

	                BL_ASSERT(pld.m_lev >= 0 && pld.m_lev <= m_gdb->finestLevel());
	                //
	                // Add it to the appropriate PBox at the appropriate level.
	                //
	                particles[pld.m_lev][std::make_pair(pld.m_grid, pld.m_tile)].push_back(p);
		    }
	        }
	    }
        }
    }

    if (ParallelDescriptor::IOProcessor() && m_verbose)
    {
        std::cout << "Done with equidistant placement" << '\n';
    }
    //
    // Let Redistribute() sort out where the particles belong.
    //
    Redistribute();

    if (ParallelDescriptor::IOProcessor() && m_verbose)
    {
        std::cout << "Redistribute done" << '\n';
    }

    BL_ASSERT(OK());
    //
    // FIXME: Will we ever need initial particles in grids deeper than 0?!
    //
    ParticleLevel& pmap = this->GetParticles(0);
    //
    // Make sure, that mf and m_gdb.boxArray(0) are defined on the same boxarray.
    //
    ParticleLocData pld;
     for (auto& kv : pmap) {
        const int        grid    = kv.first.first;
        AoS&             pbox    = kv.second.GetArrayOfStructs();
        const int        n       = pbox.size();
        const FArrayBox& dfab    = mf[grid];

#ifdef _OPENMP
#pragma omp parallel for private(pld)
#endif
        for (int i = 0; i < n; i++)
        {
            ParticleType& p = pbox[i];

            if (p.id() <= 0) continue;

            Real disp[BL_SPACEDIM];
            //
	    // Do CIC interpolation onto the particle positions.
	    // For CIC we need one ghost cell!
            //
            ParticleType::GetGravity(dfab, m_gdb->Geom(0), p, disp);

            D_TERM(p.pos(0) += len[0]*disp[0];,
                   p.pos(1) += len[1]*disp[1];,
                   p.pos(2) += len[2]*disp[2];);
            //
            // Note: m_data[0] is mass, 1 is v_x, ...
            //
            D_TERM(p.rdata(1) = vel_fac[0]*disp[0];,
                   p.rdata(2) = vel_fac[1]*disp[1];,
                   p.rdata(3) = vel_fac[2]*disp[2];);

            if (!this->Where(p, pld))
            {
	        this->PeriodicShift(p);

                if (!this->Where(p, pld))
                    amrex::Abort("DarkMatterParticleContainer::InitCosmo(): invalid particle");
	    }

            this->Reset(p, true);
        }
    }
    //
    // Let Redistribute() sort out where the particles now belong.
    //
    Redistribute();

    if (ParallelDescriptor::IOProcessor() && m_verbose)
    {
        std::cout << "Done with particle displacement" << '\n';
    }

    if (m_verbose > 1)
    {
        Real runtime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(runtime, IOProc);

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "InitCosmo() done time: " << runtime << '\n';
        }
    }
}

/*
  Particle deposition
*/

void
DarkMatterParticleContainer::AssignDensityAndVels (Array<std::unique_ptr<MultiFab> >& mf, int lev_min) const
{
     AssignDensity(mf, lev_min, BL_SPACEDIM+1);
}
