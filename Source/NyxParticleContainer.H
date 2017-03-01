
#ifndef _NyxParticleContainer_H_
#define _NyxParticleContainer_H_

#include "AMReX_Amr.H"
#include "AMReX_AmrLevel.H"
#include "AMReX_AmrParticles.H"

class NyxParticleContainerBase
{
public:
    virtual ~NyxParticleContainerBase() {}

    virtual void moveKickDrift (amrex::MultiFab& acceleration, int level, amrex::Real timestep, 
				amrex::Real a_old = 1.0, amrex::Real a_half = 1.0) = 0;
    virtual void moveKick      (amrex::MultiFab& acceleration, int level, amrex::Real timestep, 
				amrex::Real a_new = 1.0, amrex::Real a_half = 1.0,
				int start_comp_for_accel = -1) = 0;
    virtual void Redistribute (int  lev_min              = 0, 
        		       int  nGrow                = 0) = 0;
    virtual void RemoveParticlesAtLevel (int level) = 0;
    virtual amrex::Real sumParticleMass (int level) const = 0;
    virtual void AssignDensitySingleLevel (amrex::MultiFab& mf, int level, int ncomp=1,
					   int particle_lvl_offset = 0) const = 0;
    virtual void AssignDensity (amrex::Array<std::unique_ptr<amrex::MultiFab> >& mf, int lev_min = 0, int ncomp = 1,
				int finest_level = -1) const = 0;
};

template <int NR, int NI=0, int NA=0>
class NyxParticleContainer
    : public amrex::AmrParticleContainer<NR,NI,NA>,
      public NyxParticleContainerBase
{
public:


    NyxParticleContainer (amrex::Amr* amr)
	: amrex::AmrParticleContainer<NR,NI,NA>(amr),
	  sub_cycle(amr->subCycle())
    {}

    virtual ~NyxParticleContainer () {}

    void GetParticleVelocities (amrex::Array<amrex::Real>& part_vels);
    void SetParticleVelocities (amrex::Array<amrex::Real>& part_data);

    virtual amrex::Real sumParticleMass (int level) const override {
	return amrex::AmrParticleContainer<NR,NI,NA>::sumParticleMass(0,level);
    }

    void sumParticleMomentum (int lev, amrex::Real* mom) const;

    virtual void AssignDensitySingleLevel (amrex::MultiFab& mf, int level, int ncomp=1, int particle_lvl_offset = 0) const override
    { 
	amrex::AmrParticleContainer<NR,NI,NA>::AssignDensitySingleLevel(0, mf, level, ncomp, particle_lvl_offset);
    }
    virtual void AssignDensity (amrex::Array<std::unique_ptr<amrex::MultiFab> >& mf, int lev_min = 0, int ncomp = 1, int finest_level = -1) const override
    {
	amrex::AmrParticleContainer<NR,NI,NA>::AssignDensity(0, sub_cycle, mf, lev_min, ncomp, finest_level);
    }

    void MultiplyParticleMass (int lev, amrex::Real mult);

    amrex::Real estTimestep (amrex::MultiFab& acceleration,                int level, amrex::Real cfl) const;
    amrex::Real estTimestep (amrex::MultiFab& acceleration, amrex::Real a, int level, amrex::Real cfl) const;

    //
    // TODO: the methods should return a constraint on the timestep...
    //
    // The following two functions form a KICK DRIFT KICK scheme for integrating the motion of the particles in
    //   comoving coordinates -- these rely on CELL-BASED acceleration
    //
    virtual void moveKickDrift (amrex::MultiFab& acceleration, int level, amrex::Real timestep, 
				amrex::Real a_old = 1.0, amrex::Real a_half = 1.0) override;
    virtual void moveKick      (amrex::MultiFab& acceleration, int level, amrex::Real timestep, 
				amrex::Real a_new = 1.0, amrex::Real a_half = 1.0,
				int start_comp_for_accel = -1) override
    { 
	amrex::AmrParticleContainer<NR,NI,NA>::moveKick(acceleration, level, timestep,
					      a_new, a_half, start_comp_for_accel);
    }

    virtual void Redistribute (int  lev_min              = 0, 
        		       int  nGrow                = 0) override
    {
        amrex::AmrParticleContainer<NR,NI,NA>::Redistribute(lev_min, nGrow);
    }

    virtual void RemoveParticlesAtLevel (int level) override
    {
	amrex::AmrParticleContainer<NR,NI,NA>::RemoveParticlesAtLevel(level);
    }
		     
    typedef amrex::Particle<NR, NI> ParticleType;
    using AoS = typename amrex::ParticleContainer<NR,NI,NA>::AoS;
    using ParticleLevel = typename amrex::ParticleContainer<NR,NI,NA>::ParticleLevel;

protected:
    bool sub_cycle;
};

template <int NR, int NI, int NA>
void
NyxParticleContainer<NR,NI,NA>::GetParticleVelocities (amrex::Array<amrex::Real>& part_data)
{
    BL_PROFILE("NyxParticleContainer<NR,NI,NA>::GetParticleVelocities()");
    // This assumes that the mass/charge is stored in the first position 
    //      in the particle data, followed by the velocity components
    int start_comp = 1;
    int   num_comp = BL_SPACEDIM;
    this->GetParticleData(part_data,1,BL_SPACEDIM);
}

template <int NR, int NI, int NA>
void
NyxParticleContainer<NR,NI,NA>::SetParticleVelocities (amrex::Array<amrex::Real>& part_data)
{
    BL_PROFILE("NyxParticleContainer<NR,NI,NA>::SetParticleVelocities()");
   // This gives us the starting point into the part_data array
   // If only one processor (or no MPI), then that's all we need
   int cnt = 0;

#if BL_USE_MPI
   amrex::Array<long> cnts(amrex::ParallelDescriptor::NProcs());

   // This returns the number of particles on this processor
   long lcnt = this->TotalNumberOfParticles(true,true);

   // This accumulates the "lcnt" values into "cnts"
   MPI_Gather(&lcnt,1,              
              amrex::ParallelDescriptor::Mpi_typemap<long>::type(),
              cnts.dataPtr(),
              1,
              amrex::ParallelDescriptor::Mpi_typemap<long>::type(),
              amrex::ParallelDescriptor::IOProcessorNumber(),
              amrex::ParallelDescriptor::Communicator());

   amrex::ParallelDescriptor::Bcast(cnts.dataPtr(), cnts.size(), amrex::ParallelDescriptor::IOProcessorNumber());

   for (int iproc = 0; iproc < amrex::ParallelDescriptor::MyProc(); iproc++)
       cnt += cnts[iproc];

   // Each particle takes up (BL_SPACEDIM) Reals
   cnt*= (BL_SPACEDIM);
#endif

   // This is the total number of particles on *all* processors
   long npart = this->TotalNumberOfParticles(true,false);

   // Velocities
   if (part_data.size() != npart*(BL_SPACEDIM))
       amrex::Abort("Sending in wrong size part_data to SetParticleVelocities");

   for (int lev = 0; lev <= this->m_gdb->finestLevel(); lev++)
   {
       ParticleLevel& pmap = this->GetParticles(lev);

       for (typename ParticleLevel::iterator pmap_it = pmap.begin(), pmapEnd = pmap.end(); pmap_it != pmapEnd; ++pmap_it)
       {
           AoS&     pbx = pmap_it->second.GetAoS();
           const int n    = pbx.size();
    
           for (int i = 0; i < n; i++)
           {
              ParticleType& p = pbx[i];
              if (p.id() > 0)
              {
                  // Load velocities
                  for (int d=0; d < BL_SPACEDIM; d++)
                      p.rdata(d+1) = part_data[cnt+d];

                  // Update counter
                  cnt += BL_SPACEDIM;
              }
           }
       }
    }
}


//
// Assumes mass is in rdata(0), vx in rdata(1), ...!
// dim defines the cartesian direction in which the momentum is summed, x is 0, y is 1, ...
//

template <int NR, int NI, int NA>
void
NyxParticleContainer<NR,NI,NA>::sumParticleMomentum (int   lev,
					       amrex::Real* mom) const
{
    BL_PROFILE("NyxParticleContainer<NR,NI,NA>::sumParticleMomentum()");
    BL_ASSERT(NR >= BL_SPACEDIM+1);
    BL_ASSERT(lev >= 0 && lev < this->GetParticles().size());

    const ParticleLevel& pmap = this->GetParticles(lev);

    D_TERM(mom[0] = 0;, mom[1] = 0;, mom[2] = 0;);

    for (typename ParticleLevel::const_iterator pmap_it = pmap.begin(), pmapEnd = pmap.end(); pmap_it != pmapEnd; ++pmap_it)
    {
        const AoS& pbox = pmap_it->second.GetAoS();
        const int   n    = pbox.size();

        amrex::Real mom_0 = 0, mom_1 = 0, mom_2 = 0;

#ifdef _OPENMP
#pragma omp parallel for reduction(+:mom_0,mom_1,mom_2)
#endif
        for (int i = 0; i < n; i++)
        {
            const ParticleType& p = pbox[i];

            if (p.id() > 0)
            {
                D_TERM(mom_0 += p.rdata(0) * p.rdata(1);,
                       mom_1 += p.rdata(0) * p.rdata(2);,
                       mom_2 += p.rdata(0) * p.rdata(3););
            }
        }
        
        D_TERM(mom[0] += mom_0;, mom[1] += mom_1;, mom[2] += mom_2;);
    }

    amrex::ParallelDescriptor::ReduceRealSum(mom,BL_SPACEDIM);
}


template <int NR, int NI, int NA>
amrex::Real
NyxParticleContainer<NR,NI,NA>::estTimestep (amrex::MultiFab&       acceleration,
				       int             lev,
				       amrex::Real            cfl) const
{
    return estTimestep(acceleration,1.0,lev,cfl);
}

template <int NR, int NI, int NA>
amrex::Real
NyxParticleContainer<NR,NI,NA>::estTimestep (amrex::MultiFab&       acceleration,
				       amrex::Real            a,
				       int             lev,
				       amrex::Real            cfl) const
{
    BL_PROFILE("NyxParticleContainer<NR,NI,NA>::estTimestep(lev)");
    amrex::Real            dt               = 1e50;
    BL_ASSERT(NR >= BL_SPACEDIM+1);
    BL_ASSERT(lev >= 0);

    if (this->GetParticles().size() == 0)
        return dt;

    const amrex::Real      strttime         = amrex::ParallelDescriptor::second();
    const amrex::Geometry& geom             = this->m_gdb->Geom(lev);
    const amrex::Real*     dx               = geom.CellSize();
    const amrex::Real      adx[BL_SPACEDIM] = { D_DECL(a*dx[0],a*dx[1],a*dx[2]) };
    const ParticleLevel&   pmap             = this->GetParticles(lev);
    int             tnum             = 1;

#ifdef _OPENMP
    tnum = omp_get_max_threads();
#endif

    amrex::Array<amrex::Real> ldt(tnum,1e50);

    long num_particles_at_level = 0;

    amrex::MultiFab* ac_pointer;
    if (this->OnSameGrids(lev, acceleration))
    {
        ac_pointer = 0;
    }
    else 
    {
        ac_pointer = new amrex::MultiFab(this->m_gdb->ParticleBoxArray(lev), 
					 this->m_gdb->ParticleDistributionMap(lev),
					 acceleration.nComp(), acceleration.nGrow());
				  
        ac_pointer->copy(acceleration,0,0,BL_SPACEDIM);
        ac_pointer->FillBoundary(); // DO WE NEED GHOST CELLS FILLED ???
    }

    for (typename ParticleLevel::const_iterator pmap_it = pmap.begin(), pmapEnd = pmap.end(); pmap_it != pmapEnd; ++pmap_it)
    {
        const int        grid = pmap_it->first.first;
        const AoS&       pbox = pmap_it->second.GetAoS();
        const int        n    = pbox.size();
        const amrex::FArrayBox& gfab = (ac_pointer) ? (*ac_pointer)[grid] : acceleration[grid];

        num_particles_at_level += n;

#ifdef _OPENMP
#pragma omp parallel for
#endif

        for (int i = 0; i < n; i++)
        {
            const ParticleType& p = pbox[i];

            if (p.id() <= 0) continue;

            amrex::IntVect cell = this->Index(p, lev);

            const amrex::Real mag_vel_over_dx[BL_SPACEDIM] = { D_DECL(std::abs(p.rdata(1))/adx[0],
                                                                      std::abs(p.rdata(2))/adx[1],
                                                                      std::abs(p.rdata(3))/adx[2]) };

            amrex::Real max_mag_vel_over_dx = mag_vel_over_dx[0];

#if (BL_SPACEDIM > 1)
            max_mag_vel_over_dx = std::max(mag_vel_over_dx[1], max_mag_vel_over_dx);
#endif
#if (BL_SPACEDIM > 2)
            max_mag_vel_over_dx = std::max(mag_vel_over_dx[2], max_mag_vel_over_dx);
#endif
            amrex::Real dt_part = (max_mag_vel_over_dx > 0) ? (cfl / max_mag_vel_over_dx) : 1e50;

            const amrex::Real aval[BL_SPACEDIM] = { D_DECL(gfab(cell,0),
                                                           gfab(cell,1),
                                                           gfab(cell,2)) };

            const amrex::Real mag_accel = sqrt(D_TERM(aval[0]*aval[0],
                                                      + aval[1]*aval[1],
                                                      + aval[2]*aval[2]));
            if (mag_accel > 0)
                dt_part = std::min( dt_part, 1/std::sqrt(mag_accel/dx[0]) );

            int tid = 0;

#ifdef _OPENMP
            tid = omp_get_thread_num();
#endif
            ldt[tid] = std::min(dt_part, ldt[tid]);
        }
    }

    if (ac_pointer) delete ac_pointer;

    for (int i = 0; i < ldt.size(); i++)
        dt = std::min(dt, ldt[i]);

    amrex::ParallelDescriptor::ReduceRealMin(dt);
    //
    // Set dt negative if there are no particles at this level.
    //
    amrex::ParallelDescriptor::ReduceLongSum(num_particles_at_level);

    if (num_particles_at_level == 0) dt = -1.e50;

    if (this->m_verbose > 1)
    {
        amrex::Real stoptime = amrex::ParallelDescriptor::second() - strttime;

        amrex::ParallelDescriptor::ReduceRealMax(stoptime,amrex::ParallelDescriptor::IOProcessorNumber());

        if (amrex::ParallelDescriptor::IOProcessor())
        {
            std::cout << "NyxParticleContainer<NR,NI,NA>::estTimestep() time: " << stoptime << '\n';
        }
    }

    return dt;
}


template <int NR, int NI, int NA>
void
NyxParticleContainer<NR,NI,NA>::MultiplyParticleMass (int lev, amrex::Real mult)
{
    BL_PROFILE("NyxParticleContainer<NR,NI,NA>::MultiplyParticleMass()");
   BL_ASSERT(lev == 0);

   ParticleLevel& pmap = this->GetParticles(lev);

   for (typename ParticleLevel::iterator pmap_it = pmap.begin(), pmapEnd = pmap.end(); pmap_it != pmapEnd; ++pmap_it)
   {
       AoS&       pbx = pmap_it->second.GetAoS();
       const int n    = pbx.size();

#ifdef _OPENMP
#pragma omp parallel for
#endif
       for (int i = 0; i < n; i++)
       {
          ParticleType& p = pbx[i];
          if (p.id() > 0)
          {
              //
              // Note: rdata(0) is mass, ...
              //
              p.rdata(0) *= mult;
          }
       }
   }
}

//
// This version takes as input the acceleration vector at cell centers
//
template <int NR, int NI, int NA>
void
NyxParticleContainer<NR,NI,NA>::moveKickDrift (amrex::MultiFab&       acceleration,
					 int             lev,
					 amrex::Real            dt,
					 amrex::Real            a_old,
					 amrex::Real            a_half) 
{
    BL_PROFILE("ParticleContainer::moveKickDrift()");
    BL_ASSERT(NR >= BL_SPACEDIM+1);
    BL_ASSERT(lev >= 0);
    BL_ASSERT(acceleration.nGrow() >= 2);

    //If there are no particles at this level
    if (lev >= this->GetParticles().size())
        return;

    const amrex::Real strttime      = amrex::ParallelDescriptor::second();
    const amrex::Real half_dt       = amrex::Real(0.5) * dt;
    const amrex::Real a_half_inv    = 1 / a_half;
    const amrex::Real dt_a_half_inv = dt * a_half_inv;
    ParticleLevel&    pmap          = this->GetParticles(lev);

    amrex::MultiFab* ac_pointer;
    if (this->OnSameGrids(lev, acceleration))
    {
        ac_pointer = &acceleration;
    }
    else
    {
        ac_pointer = new amrex::MultiFab(this->m_gdb->ParticleBoxArray(lev),
					 this->m_gdb->ParticleDistributionMap(lev),
					 acceleration.nComp(),acceleration.nGrow());
        for (amrex::MFIter mfi(*ac_pointer); mfi.isValid(); ++mfi)
            ac_pointer->setVal(0.);
        ac_pointer->copy(acceleration,0,0,acceleration.nComp());
        ac_pointer->FillBoundary(); // DO WE NEED GHOST CELLS FILLED ???
    }

    for (auto& kv : pmap) {
        const int        grid = kv.first.first;
        AoS&             pbox = kv.second.GetAoS();
        const int        n    = pbox.size();
        const amrex::FArrayBox& gfab = (*ac_pointer)[grid];

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < n; i++)
        {
            ParticleType& p = pbox[i];

            if (p.id() <= 0) continue;

            //
            // note: m_data[0] is mass, 1 is v_x, ...
            //
            amrex::Real grav[BL_SPACEDIM];

            ParticleType::GetGravity(gfab, this->m_gdb->Geom(lev), p, grav);
            //
            // First update (a u)^half = (a u)^old + dt/2 grav^old
            //
            D_TERM(p.rdata(1) *= a_old;,
                   p.rdata(2) *= a_old;,
                   p.rdata(3) *= a_old;);

            //
            // Add adot/a and gravitational updates.
            //
            D_TERM(p.rdata(1) += half_dt * grav[0];,
                   p.rdata(2) += half_dt * grav[1];,
                   p.rdata(3) += half_dt * grav[2];);

            D_TERM(p.rdata(1) *= a_half_inv;,
                   p.rdata(2) *= a_half_inv;,
                   p.rdata(3) *= a_half_inv;);
            //
            // Now update x^new = x^old + dt grav^half / a^half
            //

            D_TERM(p.pos(0) += dt_a_half_inv * p.rdata(1);,
                   p.pos(1) += dt_a_half_inv * p.rdata(2);,
                   p.pos(2) += dt_a_half_inv * p.rdata(3););
        }
    }

    if (ac_pointer != &acceleration) delete ac_pointer;
    
    if (lev > 0 && sub_cycle)
    {
        for (auto& kv : pmap) {
            AoS&  pbox       = kv.second.GetAoS();
            const int   n    = pbox.size();

            amrex::ParticleLocData pld; 
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (int i = 0; i < n; i++)
            {
                ParticleType& p = pbox[i];
                if (p.id() <= 0) continue;

                // Move the particle to the proper ghost cell. 
                //      and remove any *ghost* particles that have gone too far
                // Note that this should only negate ghost particles, not real particles.
                if (!this->SingleLevelGrownWhere(p, pld, lev, acceleration.nGrow()-2))
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

    if (this->m_verbose > 1)
    {
        amrex::Real stoptime = amrex::ParallelDescriptor::second() - strttime;

        amrex::ParallelDescriptor::ReduceRealMax(stoptime,amrex::ParallelDescriptor::IOProcessorNumber());

        if (amrex::ParallelDescriptor::IOProcessor())
        {
            std::cout << "NyxParticleContainer<NR,NI,NA>::moveKickDrift() time: " << stoptime << '\n';
        }
    }
}


// We make DarkMatterParticleContainer a class instead of a typedef so that
//    we can have Nyx-specific functions here instead of in BoxLib
class DarkMatterParticleContainer
    : public NyxParticleContainer<1+BL_SPACEDIM>
{
public:
    DarkMatterParticleContainer (amrex::Amr* amr)
        : NyxParticleContainer<1+BL_SPACEDIM>(amr)
    { }

    virtual ~DarkMatterParticleContainer () {}

    void InitCosmo  (amrex::MultiFab& mf, const amrex::Real vel_fac[], const amrex::Array<int> n_part, const amrex::Real particleMass);
    void InitCosmo  (amrex::MultiFab& mf, const amrex::Real vel_fac[], const amrex::Array<int> n_part,
                     const amrex::Real particleMass, const amrex::Real shift[]);
    void InitCosmo1ppc (amrex::MultiFab& mf, const amrex::Real vel_fac[], const amrex::Real particleMass);
    void InitCosmo1ppcMultiLevel(amrex::MultiFab& mf, const amrex::Real disp_fac[], const amrex::Real vel_fac[], 
                                 const amrex::Real particleMass, int disp_idx, int vel_idx, 
                                 amrex::BoxArray &baWhereNot, int lev, int nlevs);

    void AssignDensityAndVels (amrex::Array<std::unique_ptr<amrex::MultiFab> >& mf, int lev_min = 0) const;

};

// We make NeutrinoParticleContainer a class instead of a typedef so that
//    we can have Nyx-specific functions here instead of in BoxLib
class NeutrinoParticleContainer
    : public NyxParticleContainer<2+BL_SPACEDIM>
{

    int         m_relativistic; // if 1 then we weight the mass by gamma in AssignDensity*
    amrex::Real        m_csq;          // the square of the speed of light -- used to compute relativistic effects

public:
    NeutrinoParticleContainer (amrex::Amr* amr)
        : NyxParticleContainer<2+BL_SPACEDIM>(amr)
    { }

    virtual ~NeutrinoParticleContainer () {}

    void SetRelativistic (int relativistic) { m_relativistic = relativistic; }

    void SetCSquared (amrex::Real csq) { m_csq = csq; }

    void AssignDensity (amrex::Array<std::unique_ptr<amrex::MultiFab> >& mf, int lev_min = 0, int ncomp = 1, int finest_level = -1) const;

    void AssignRelativisticDensitySingleLevel (amrex::MultiFab& mf, int level, int ncomp=1, int particle_lvl_offset = 0) const;
    
    void AssignRelativisticDensity (amrex::Array<std::unique_ptr<amrex::MultiFab> >& mf, int lev_min = 0, int ncomp = 1, int finest_level = -1) const;

};

#endif /*_NyxParticleContainer_H_*/
