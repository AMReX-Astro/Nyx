#include "NyxParticleContainer.H"

void
DarkMatterParticleContainer::InitCosmo1ppcMultiLevel(
                        MultiFab& mf, const Real disp_fac[], const Real vel_fac[], 
                        const Real particleMass, int disp_idx, int vel_idx, 
                        BoxArray &baWhereNot, int lev)
{
    BL_PROFILE("ParticleContainer<N>::InitCosmo1ppcMultiLevel()");
    const int       MyProc   = ParallelDescriptor::MyProc();
    const Geometry& geom     = m_gdb->Geom(lev);
    const Real*     dx       = geom.CellSize();

    static Array<int> calls;

    calls.resize(m_gdb->initialBaLevels()+1);

    calls[lev]++;

    if (calls[lev] > 1) return;

    m_particles.reserve(15);  // So we don't ever have to do any copying on a resize.

    m_particles.resize(m_gdb->initialBaLevels()+1);

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
        const int  *fab_lo = mfi.validbox().loVect();
        const int  *fab_hi = mfi.validbox().hiVect();

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
	                p.m_pos[n] = geom.ProbLo(n) + 
                                      (indices[n]+Real(0.5))*dx[n];
			if(disp[n]*disp_fac[n]>dx[n]/2.0)
			  outcount[n]++;
			if(disp[n]*disp_fac[n]<-dx[n]/2.0)
			  outcountminus[n]++;
			mean_disp[n]+=fabs(disp[n]);
			//
                        // then add the displacement (input values weighted by domain length).
                        //
	                p.m_pos[n] += disp[n] * disp_fac[n];

                        //
			// Set the velocities.
                        //
                        vel[n] = myFab(indices,vel_idx+n);
	                p.m_data[n+1] = vel[n] * vel_fac[n];
	            }
                    //
		    // Set the mass of the particle from the input value.
                    //
	            p.m_data[0] = particleMass;
	            p.m_id      = ParticleBase::NextID();
	            p.m_cpu     = MyProc;
	
	            if (!ParticleBase::Where(p,m_gdb))
                    {
      		        ParticleBase::PeriodicShift(p,m_gdb);

                        if (!ParticleBase::Where(p,m_gdb))
                            BoxLib::Abort("ParticleContainer<N>::InitCosmo1ppcMultiLevel():invalid particle");
                    }

		    BL_ASSERT(p.m_lev >= 0 && p.m_lev <= m_gdb->finestLevel());
		    //handle particles that ran out of this level into a finer one. 
		    if (baWhereNot.contains(p.m_cell))
		    {
		      outside_counter++;
		      ParticleType newp[8];
		      for (int i=0;i<8;i++)
		      {
			newp[i].m_data[0] = particleMass/8.0;
			newp[i].m_id = ParticleBase::NextID();
			newp[i].m_cpu     = MyProc;
			for (int dim=0;dim<BL_SPACEDIM;dim++)
			{
			  newp[i].m_pos[dim]=p.m_pos[dim]+(2*((i/(1 << dim)) % 2)-1)*dx[dim]/4.0;
			  newp[i].m_data[dim+1]=p.m_data[dim+1];
			  
			  
			}
		      
			if (!ParticleBase::Where(newp[i],m_gdb))
			{
			    ParticleBase::PeriodicShift(newp[i],m_gdb);
    
			    if (!ParticleBase::Where(newp[i],m_gdb))
                                BoxLib::Abort("ParticleContainer<N>::InitCosmo1ppcMultiLevel():invalid particle");
			}
			m_particles[newp[i].m_lev][newp[i].m_grid].push_back(newp[i]);
		      }
		      
		    }
	            
	            //
	            // Add it to the appropriate PBox at the appropriate level.
	            //
		    else
		      m_particles[p.m_lev][p.m_grid].push_back(p);
                }
            }
        }
    }
    Redistribute();
}

void
DarkMatterParticleContainer::InitCosmo1ppc(MultiFab& mf, const Real vel_fac[], const Real particleMass)
{
    BL_PROFILE("ParticleContainer<N>::InitCosmo1ppc()");
    const int       MyProc   = ParallelDescriptor::MyProc();
    const Geometry& geom     = m_gdb->Geom(0);
    const Real*     dx       = geom.CellSize();

    m_particles.reserve(15);  // So we don't ever have to do any copying on a resize.

    m_particles.resize(m_gdb->finestLevel()+1);

    for (int lev = 0; lev < m_particles.size(); lev++)
    {
        BL_ASSERT(m_particles[lev].empty());
    }

    ParticleType p;
    Real         disp[BL_SPACEDIM];
    const Real   len[BL_SPACEDIM] = { D_DECL(geom.ProbLength(0),
                                             geom.ProbLength(1),
                                             geom.ProbLength(2)) };
    //
    // The grid should be initialized according to the ics...
    //
    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        FArrayBox&  myFab  = mf[mfi];
        const int  *fab_lo = mfi.validbox().loVect();
        const int  *fab_hi = mfi.validbox().hiVect();

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
	                p.m_pos[n] = geom.ProbLo(n) + 
                                      (indices[n]+Real(0.5))*dx[n] +
		                      disp[n] * len[n];
                        //
			// Set the velocities.
                        //
	                p.m_data[n+1] = disp[n] * vel_fac[n];
	            }
                    //
		    // Set the mass of the particle from the input value.
                    //
	            p.m_data[0] = particleMass;
	            p.m_id      = ParticleBase::NextID();
	            p.m_cpu     = MyProc;
	
	            if (!ParticleBase::Where(p,m_gdb))
                    {
      		        ParticleBase::PeriodicShift(p,m_gdb);

                        if (!ParticleBase::Where(p,m_gdb))
                            BoxLib::Abort("ParticleContainer<N>::InitCosmo1ppc(): invalid particle");
		    }

	            BL_ASSERT(p.m_lev >= 0 && p.m_lev <= m_gdb->finestLevel());
	            //
	            // Add it to the appropriate PBox at the appropriate level.
	            //
	            m_particles[p.m_lev][p.m_grid].push_back(p);
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
    BL_PROFILE("ParticleContainer<N>::InitCosmo()");
    const int       MyProc   = ParallelDescriptor::MyProc();
    const int       IOProc   = ParallelDescriptor::IOProcessorNumber();
    const Real      strttime = ParallelDescriptor::second();
    const Geometry& geom     = m_gdb->Geom(0);

    m_particles.reserve(15);  // So we don't ever have to do any copying on a resize.

    m_particles.resize(m_gdb->finestLevel()+1);

    for (int lev = 0; lev < m_particles.size(); lev++)
    {
        BL_ASSERT(m_particles[lev].empty());
    }

    const Real len[BL_SPACEDIM] = { D_DECL(geom.ProbLength(0),
                                           geom.ProbLength(1),
                                           geom.ProbLength(2)) };
    //
    // Print the grids as a sanity check.
    //
    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        const int  *fab_lo = mfi.validbox().loVect();
        const int  *fab_hi = mfi.validbox().hiVect();
        if (mfi.validbox().isEmpty())
        {
           std::cout << "...bad grid lo " << fab_lo[0] << " " << fab_lo[1] << " " << fab_lo[2] << '\n';
           std::cout << "...bad grid hi " << fab_hi[0] << " " << fab_hi[1] << " " << fab_hi[2] << '\n';
           BoxLib::Error("Empty box in InitCosmo ");
        }
        if (!geom.Domain().contains(mfi.validbox()))
        {
           std::cout << "...bad grid lo " << fab_lo[0] << " " << fab_lo[1] << " " << fab_lo[2] << '\n';
           std::cout << "...bad grid hi " << fab_hi[0] << " " << fab_hi[1] << " " << fab_hi[2] << '\n';
           BoxLib::Error("Box in InitCosmo not contained in domain");
        }
    }
    //
    // We will need one ghost cell, so check wether we have one.
    //
    if (mf.nGrow() < 1)
        BoxLib::Abort("ParticleContainer<N>::InitCosmo: mf needs at least one correctly filled ghost zone!");

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
                        D_TERM(p.m_pos[0] = pos[0];,
                               p.m_pos[1] = pos[1];,
                               p.m_pos[2] = pos[2];);
                        //
		        // Set the mass of the particle.
                        //
	                p.m_data[0] = particleMass;
	                p.m_id      = ParticleBase::NextID();
	                p.m_cpu     = MyProc;

	                if (!ParticleBase::Where(p,m_gdb))
                        {
      		            ParticleBase::PeriodicShift(p,m_gdb);

                            if (!ParticleBase::Where(p,m_gdb))
                                BoxLib::Abort("ParticleContainer<N>::InitCosmo(): invalid particle");
		        }

	                BL_ASSERT(p.m_lev >= 0 && p.m_lev <= m_gdb->finestLevel());
	                //
	                // Add it to the appropriate PBox at the appropriate level.
	                //
	                m_particles[p.m_lev][p.m_grid].push_back(p);
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
    Redistribute(true);

    if (ParallelDescriptor::IOProcessor() && m_verbose)
    {
        std::cout << "Redistribute done" << '\n';
    }

    BL_ASSERT(OK());
    BL_ASSERT(N >= BL_SPACEDIM+1);
    //
    // FIXME: Will we ever need initial particles in grids deeper than 0?!
    //
    PMap& pmap = m_particles[0];
    //
    // Make sure, that mf and m_gdb.boxArray(0) are defined on the same boxarray.
    //
    for (typename PMap::iterator pmap_it = pmap.begin(), pmapEnd = pmap.end(); pmap_it != pmapEnd; ++pmap_it)
    {
        const int        grid    = pmap_it->first;
        PBox&            pbox    = pmap_it->second;
        const int        n       = pbox.size();
        const FArrayBox& dfab    = mf[grid];

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < n; i++)
        {
            ParticleType& p = pbox[i];

            if (p.m_id <= 0) continue;

            BL_ASSERT(p.m_grid == grid);

            Real disp[BL_SPACEDIM];
            //
	    // Do CIC interpolation onto the particle positions.
	    // For CIC we need one ghost cell!
            //
            ParticleBase::GetGravity(dfab,m_gdb->Geom(p.m_lev),p,disp);

            D_TERM(p.m_pos[0] += len[0]*disp[0];,
                   p.m_pos[1] += len[1]*disp[1];,
                   p.m_pos[2] += len[2]*disp[2];);
            //
            // Note: m_data[0] is mass, 1 is v_x, ...
            //
            D_TERM(p.m_data[1] = vel_fac[0]*disp[0];,
                   p.m_data[2] = vel_fac[1]*disp[1];,
                   p.m_data[3] = vel_fac[2]*disp[2];);

            if (!ParticleBase::Where(p,m_gdb))
            {
	        ParticleBase::PeriodicShift(p,m_gdb);

                if (!ParticleBase::Where(p,m_gdb))
                    BoxLib::Abort("ParticleContainer<N>::InitCosmo(): invalid particle");
	    }

            ParticleBase::Reset(p,m_gdb,true);
        }
    }
    //
    // Let Redistribute() sort out where the particles now belong.
    //
    Redistribute(true);

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
