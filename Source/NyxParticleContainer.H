
#ifndef _NyxParticleContainer_H_
#define _NyxParticleContainer_H_

#include "AMReX_Amr.H"
#include "AMReX_AmrLevel.H"
#include "AMReX_AmrParticles.H"

class NyxParticleContainerBase
{
public:
    virtual ~NyxParticleContainerBase() {}

    virtual void moveKickDrift (amrex::MultiFab& acceleration, int level, Real timestep, 
				Real a_old = 1.0, Real a_half = 1.0) = 0;
    virtual void moveKick      (amrex::MultiFab& acceleration, int level, Real timestep, 
				Real a_new = 1.0, Real a_half = 1.0,
				int start_comp_for_accel = -1) = 0;
    virtual void Redistribute (bool where_already_called = false,
			       bool full_where           = false,
			       int  lev_min              = 0, 
			       int  nGrow                = 0) = 0;
    virtual void RemoveParticlesAtLevel (int level) = 0;
    virtual Real sumParticleMass (int level) const = 0;
    virtual void AssignDensitySingleLevel (amrex::MultiFab& mf, int level, int ncomp=1,
					   int particle_lvl_offset = 0) const = 0;
    virtual void AssignDensity (amrex::Array<std::unique_ptr<amrex::MultiFab> >& mf, int lev_min = 0, int ncomp = 1,
				int finest_level = -1) const = 0;
};

template <int NR, int NI=0>
class NyxParticleContainer
    : public amrex::AmrParticleContainer<NR,NI>,
      public NyxParticleContainerBase
{
public:


    NyxParticleContainer (amrex::Amr* amr)
	: amrex::AmrParticleContainer<NR,NI>(amr),
	  sub_cycle(amr->subCycle())
    {}

    virtual ~NyxParticleContainer () {}

    void GetParticleVelocities (amrex::Array<Real>& part_vels);
    void SetParticleVelocities (amrex::Array<Real>& part_data);

    virtual Real sumParticleMass (int level) const override {
	return amrex::AmrParticleContainer<NR,NI>::sumParticleMass(0,level);
    }

    void sumParticleMomentum (int lev, Real* mom) const;

    virtual void AssignDensitySingleLevel (amrex::MultiFab& mf, int level, int ncomp=1, int particle_lvl_offset = 0) const override
    { 
	amrex::AmrParticleContainer<NR,NI>::AssignDensitySingleLevel(0, mf, level, ncomp, particle_lvl_offset);
    }
    virtual void AssignDensity (amrex::Array<std::unique_ptr<amrex::MultiFab> >& mf, int lev_min = 0, int ncomp = 1, int finest_level = -1) const override
    {
	amrex::AmrParticleContainer<NR,NI>::AssignDensity(0, sub_cycle, mf, lev_min, ncomp, finest_level);
    }

    void MultiplyParticleMass (int lev, Real mult);

    Real estTimestep (amrex::MultiFab& acceleration,         int level, Real cfl) const;
    Real estTimestep (amrex::MultiFab& acceleration, Real a, int level, Real cfl) const;

    //
    // TODO: the methods should return a constraint on the timestep...
    //
    // The following two functions form a KICK DRIFT KICK scheme for integrating the motion of the particles in
    //   comoving coordinates -- these rely on CELL-BASED acceleration
    //
    virtual void moveKickDrift (amrex::MultiFab& acceleration, int level, Real timestep, 
				Real a_old = 1.0, Real a_half = 1.0) override;
    virtual void moveKick      (amrex::MultiFab& acceleration, int level, Real timestep, 
				Real a_new = 1.0, Real a_half = 1.0,
				int start_comp_for_accel = -1) override
    { 
	amrex::AmrParticleContainer<NR,NI>::moveKick(acceleration, level, timestep,
					      a_new, a_half, start_comp_for_accel);
    }

    virtual void Redistribute (bool where_already_called = false,
			       bool full_where           = false,
			       int  lev_min              = 0, 
			       int  nGrow                = 0) override
    {
	amrex::AmrParticleContainer<NR,NI>::Redistribute(where_already_called,
						  full_where, lev_min, nGrow);
    }

    virtual void RemoveParticlesAtLevel (int level) override
    {
	amrex::AmrParticleContainer<NR,NI>::RemoveParticlesAtLevel(level);
    }
		     
    typedef amrex::Particle<NR,NI> ParticleType;
    //    typedef typename std::deque<ParticleType> PBox;
    //    typedef typename std::map<int,PBox> PMap;
    typedef typename amrex::AmrParticleContainer<NR,NI>::PBox PBox;
    typedef typename amrex::AmrParticleContainer<NR,NI>::PMap PMap;

protected:
    bool sub_cycle;
};

template <int NR, int NI>
void
NyxParticleContainer<NR,NI>::GetParticleVelocities (amrex::Array<Real>& part_data)
{
    BL_PROFILE("NyxParticleContainer<NR,NI>::GetParticleVelocities()");
    // This assumes that the mass/charge is stored in the first position 
    //      in the particle data, followed by the velocity components
    int start_comp = 1;
    int   num_comp = BL_SPACEDIM;
    this->GetParticleData(part_data,1,BL_SPACEDIM);
}

template <int NR, int NI>
void
NyxParticleContainer<NR,NI>::SetParticleVelocities (amrex::Array<Real>& part_data)
{
    BL_PROFILE("NyxParticleContainer<NR,NI>::SetParticleVelocities()");
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
       PMap& pmap = this->m_particles[lev];

       for (typename PMap::iterator pmap_it = pmap.begin(), pmapEnd = pmap.end(); pmap_it != pmapEnd; ++pmap_it)
       {
           PBox&     pbx = pmap_it->second;
           const int n    = pbx.size();
    
           for (int i = 0; i < n; i++)
           {
              ParticleType& p = pbx[i];
              if (p.m_id > 0)
              {
                  // Load velocities
                  for (int d=0; d < BL_SPACEDIM; d++)
                     p.m_data[d+1] = part_data[cnt+d];

                  // Update counter
                  cnt += BL_SPACEDIM;
              }
           }
       }
    }
}


//
// Assumes mass is in m_data[0], vx in m_dat[1], ...!
// dim defines the cartesian direction in which the momentum is summed, x is 0, y is 1, ...
//

template <int NR, int NI>
void
NyxParticleContainer<NR,NI>::sumParticleMomentum (int   lev,
					       Real* mom) const
{
    BL_PROFILE("NyxParticleContainer<NR,NI>::sumParticleMomentum()");
    BL_ASSERT(NR >= BL_SPACEDIM+1);
    BL_ASSERT(lev >= 0 && lev < this->m_particles.size());

    const PMap& pmap = this->m_particles[lev];

    D_TERM(mom[0] = 0;, mom[1] = 0;, mom[2] = 0;);

    for (typename PMap::const_iterator pmap_it = pmap.begin(), pmapEnd = pmap.end(); pmap_it != pmapEnd; ++pmap_it)
    {
        const PBox& pbox = pmap_it->second;
        const int   n    = pbox.size();

        Real mom_0 = 0, mom_1 = 0, mom_2 = 0;

#ifdef _OPENMP
#pragma omp parallel for reduction(+:mom_0,mom_1,mom_2)
#endif
        for (int i = 0; i < n; i++)
        {
            const ParticleType& p = pbox[i];

            if (p.m_id > 0)
            {
                D_TERM(mom_0 += p.m_data[0] * p.m_data[1];,
                       mom_1 += p.m_data[0] * p.m_data[2];,
                       mom_2 += p.m_data[0] * p.m_data[3];);
            }
        }
        
        D_TERM(mom[0] += mom_0;, mom[1] += mom_1;, mom[2] += mom_2;);
    }

    amrex::ParallelDescriptor::ReduceRealSum(mom,BL_SPACEDIM);
}


template <int NR, int NI>
Real
NyxParticleContainer<NR,NI>::estTimestep (amrex::MultiFab&       acceleration,
				       int             lev,
				       Real            cfl) const
{
    return estTimestep(acceleration,1.0,lev,cfl);
}

template <int NR, int NI>
Real
NyxParticleContainer<NR,NI>::estTimestep (amrex::MultiFab&       acceleration,
				       Real            a,
				       int             lev,
				       Real            cfl) const
{
    BL_PROFILE("NyxParticleContainer<NR,NI>::estTimestep(lev)");
    Real            dt               = 1e50;
    BL_ASSERT(NR >= BL_SPACEDIM+1);
    BL_ASSERT(lev >= 0);

    if (this->m_particles.size() == 0)
        return dt;

    const Real      strttime         = amrex::ParallelDescriptor::second();
    const amrex::Geometry& geom             = this->m_gdb->Geom(lev);
    const Real*     dx               = geom.CellSize();
    const Real      adx[BL_SPACEDIM] = { D_DECL(a*dx[0],a*dx[1],a*dx[2]) };
    const PMap&     pmap             = this->m_particles[lev];
    int             tnum             = 1;

#ifdef _OPENMP
    tnum = omp_get_max_threads();
#endif

    amrex::Array<Real> ldt(tnum,1e50);

    long num_particles_at_level = 0;

    amrex::MultiFab* ac_pointer;
    if (this->OnSameGrids(lev, acceleration))
    {
        ac_pointer = 0;
    }
    else 
    {
        ac_pointer = new amrex::MultiFab(this->m_gdb->ParticleBoxArray(lev), acceleration.nComp(), acceleration.nGrow(),
				  this->m_gdb->ParticleDistributionMap(lev), amrex::Fab_allocate);
        ac_pointer->copy(acceleration,0,0,BL_SPACEDIM);
        ac_pointer->FillBoundary(); // DO WE NEED GHOST CELLS FILLED ???
    }

    for (typename PMap::const_iterator pmap_it = pmap.begin(), pmapEnd = pmap.end(); pmap_it != pmapEnd; ++pmap_it)
    {
        const int        grid = pmap_it->first;
        const PBox&      pbox = pmap_it->second;
        const int        n    = pbox.size();
        const amrex::FArrayBox& gfab = (ac_pointer) ? (*ac_pointer)[grid] : acceleration[grid];

        num_particles_at_level += n;

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < n; i++)
        {
            const ParticleType& p = pbox[i];

            if (p.m_id <= 0) continue;

            BL_ASSERT(p.m_grid == grid);

            const Real mag_vel_over_dx[BL_SPACEDIM] = { D_DECL(std::abs(p.m_data[1])/adx[0],
                                                               std::abs(p.m_data[2])/adx[1],
                                                               std::abs(p.m_data[3])/adx[2]) };

            Real max_mag_vel_over_dx = mag_vel_over_dx[0];

#if (BL_SPACEDIM > 1)
            max_mag_vel_over_dx = std::max(mag_vel_over_dx[1], max_mag_vel_over_dx);
#endif
#if (BL_SPACEDIM > 2)
            max_mag_vel_over_dx = std::max(mag_vel_over_dx[2], max_mag_vel_over_dx);
#endif
            Real dt_part = (max_mag_vel_over_dx > 0) ? (cfl / max_mag_vel_over_dx) : 1e50;

            const Real aval[BL_SPACEDIM] = { D_DECL(gfab(p.m_cell,0),
                                                    gfab(p.m_cell,1),
                                                    gfab(p.m_cell,2)) };

            const Real mag_accel = sqrt(D_TERM(aval[0]*aval[0],
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
        Real stoptime = amrex::ParallelDescriptor::second() - strttime;

        amrex::ParallelDescriptor::ReduceRealMax(stoptime,amrex::ParallelDescriptor::IOProcessorNumber());

        if (amrex::ParallelDescriptor::IOProcessor())
        {
            std::cout << "NyxParticleContainer<NR,NI>::estTimestep() time: " << stoptime << '\n';
        }
    }

    return dt;
}


template <int NR, int NI>
void
NyxParticleContainer<NR,NI>::MultiplyParticleMass (int lev, Real mult)
{
    BL_PROFILE("NyxParticleContainer<NR,NI>::MultiplyParticleMass()");
   BL_ASSERT(lev == 0);

   PMap& pmap = this->m_particles[lev];

   for (typename PMap::iterator pmap_it = pmap.begin(), pmapEnd = pmap.end(); pmap_it != pmapEnd; ++pmap_it)
   {
       PBox&     pbx = pmap_it->second;
       const int n    = pbx.size();

#ifdef _OPENMP
#pragma omp parallel for
#endif
       for (int i = 0; i < n; i++)
       {
          ParticleType& p = pbx[i];
          if (p.m_id > 0)
          {
              //
              // Note: m_data[0] is mass, ...
              //
              p.m_data[0] *= mult;
          }
       }
   }
}

//
// This version takes as input the acceleration vector at cell centers
//
template <int NR, int NI>
void
NyxParticleContainer<NR,NI>::moveKickDrift (amrex::MultiFab&       acceleration,
					 int             lev,
					 Real            dt,
					 Real            a_old,
					 Real            a_half) 
{
    BL_PROFILE("ParticleContainer::moveKickDrift()");
    BL_ASSERT(NR >= BL_SPACEDIM+1);
    BL_ASSERT(lev >= 0);
    BL_ASSERT(acceleration.nGrow() >= 2);

    //If there are no particles at this level
    if (lev >= this->m_particles.size())
        return;

    const Real strttime      = amrex::ParallelDescriptor::second();
    const Real half_dt       = Real(0.5) * dt;
    const Real a_half_inv    = 1 / a_half;
    const Real dt_a_half_inv = dt * a_half_inv;
    PMap&      pmap          = this->m_particles[lev];

    amrex::MultiFab* ac_pointer;
    if (this->OnSameGrids(lev, acceleration))
    {
        ac_pointer = &acceleration;
    }
    else
    {
        ac_pointer = new amrex::MultiFab(this->m_gdb->ParticleBoxArray(lev),acceleration.nComp(),acceleration.nGrow(),
				  this->m_gdb->ParticleDistributionMap(lev),amrex::Fab_allocate);
        for (amrex::MFIter mfi(*ac_pointer); mfi.isValid(); ++mfi)
            ac_pointer->setVal(0.);
        ac_pointer->copy(acceleration,0,0,acceleration.nComp());
        ac_pointer->FillBoundary(); // DO WE NEED GHOST CELLS FILLED ???
    }

    for (typename PMap::iterator pmap_it = pmap.begin(), pmapEnd = pmap.end(); pmap_it != pmapEnd; ++pmap_it)
    {
        const int        grid = pmap_it->first;
        PBox&            pbox = pmap_it->second;
        const int        n    = pbox.size();
        const amrex::FArrayBox& gfab = (*ac_pointer)[grid];

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < n; i++)
        {
            ParticleType& p = pbox[i];

            if (p.m_id <= 0) continue;

            BL_ASSERT(p.m_grid == grid);
            //
            // note: m_data[0] is mass, 1 is v_x, ...
            //
            Real grav[BL_SPACEDIM];

            amrex::ParticleBase::GetGravity(gfab, this->m_gdb->Geom(p.m_lev), p, grav);
            //
            // First update (a u)^half = (a u)^old + dt/2 grav^old
            //
            D_TERM(p.m_data[1] *= a_old;,
                   p.m_data[2] *= a_old;,
                   p.m_data[3] *= a_old;);

            //
            // Add adot/a and gravitational updates.
            //
            D_TERM(p.m_data[1] += half_dt * grav[0];,
                   p.m_data[2] += half_dt * grav[1];,
                   p.m_data[3] += half_dt * grav[2];);

            D_TERM(p.m_data[1] *= a_half_inv;,
                   p.m_data[2] *= a_half_inv;,
                   p.m_data[3] *= a_half_inv;);
            //
            // Now update x^new = x^old + dt grav^half / a^half
            //

            D_TERM(p.m_pos[0] += dt_a_half_inv * p.m_data[1];,
                   p.m_pos[1] += dt_a_half_inv * p.m_data[2];,
                   p.m_pos[2] += dt_a_half_inv * p.m_data[3];);
        }
    }

    if (ac_pointer != &acceleration) delete ac_pointer;

    if (lev > 0 && sub_cycle)
    {
        for (typename PMap::iterator pmap_it = pmap.begin(), pmapEnd = pmap.end(); pmap_it != pmapEnd; ++pmap_it)
        {
            PBox&            pbox = pmap_it->second;
            const int        n    = pbox.size();
 
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (int i = 0; i < n; i++)
            {
                ParticleType& p = pbox[i];
                if (p.m_id <= 0) continue;

                // Move the particle to the proper ghost cell. 
                //      and remove any *ghost* particles that have gone too far
                // Note that this should only negate ghost particles, not real particles.
                if (!amrex::ParticleBase::RestrictedWhere(p,this->m_gdb, acceleration.nGrow()-2))
                {
                    // Assert that the particle being removed is a ghost particle;
                    // the ghost particle is no longer in relevant ghost cells for this grid.
                    if (p.m_id == amrex::GhostParticleID)
                    {
                        p.m_id = -1;
                    }
                    else
                    {
                        std::cout << "Oops -- removing particle " << p.m_id << std::endl;
                        amrex::Error("Trying to get rid of a non-ghost particle in moveKickDrift");
                    }
                }
            }
        }
    }

    if (this->m_verbose > 1)
    {
        Real stoptime = amrex::ParallelDescriptor::second() - strttime;

        amrex::ParallelDescriptor::ReduceRealMax(stoptime,amrex::ParallelDescriptor::IOProcessorNumber());

        if (amrex::ParallelDescriptor::IOProcessor())
        {
            std::cout << "NyxParticleContainer<NR,NI>::moveKickDrift() time: " << stoptime << '\n';
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

    void InitCosmo  (amrex::MultiFab& mf, const Real vel_fac[], const amrex::Array<int> n_part, const Real particleMass);
    void InitCosmo  (amrex::MultiFab& mf, const Real vel_fac[], const amrex::Array<int> n_part,
                     const Real particleMass, const Real shift[]);
    void InitCosmo1ppc (amrex::MultiFab& mf, const Real vel_fac[], const Real particleMass);
    void InitCosmo1ppcMultiLevel(amrex::MultiFab& mf, const Real disp_fac[], const Real vel_fac[], 
                                 const Real particleMass, int disp_idx, int vel_idx, 
                                 amrex::BoxArray &baWhereNot, int lev, int nlevs);

    void AssignDensityAndVels (amrex::Array<std::unique_ptr<amrex::MultiFab> >& mf, int lev_min = 0) const;

};

// We make NeutrinoParticleContainer a class instead of a typedef so that
//    we can have Nyx-specific functions here instead of in BoxLib
class NeutrinoParticleContainer
    : public NyxParticleContainer<2+BL_SPACEDIM>
{

    int         m_relativistic; // if 1 then we weight the mass by gamma in AssignDensity*
    Real        m_csq;          // the square of the speed of light -- used to compute relativistic effects

public:
    NeutrinoParticleContainer (amrex::Amr* amr)
        : NyxParticleContainer<2+BL_SPACEDIM>(amr)
    { }

    virtual ~NeutrinoParticleContainer () {}

    void SetRelativistic (int relativistic) { m_relativistic = relativistic; }

    void SetCSquared (Real csq) { m_csq = csq; }

    void AssignDensity (amrex::Array<std::unique_ptr<amrex::MultiFab> >& mf, int lev_min = 0, int ncomp = 1, int finest_level = -1) const;

    void AssignRelativisticDensitySingleLevel (amrex::MultiFab& mf, int level, int ncomp=1, int particle_lvl_offset = 0) const;
    
    void AssignRelativisticDensity (amrex::Array<std::unique_ptr<amrex::MultiFab> >& mf, int lev_min = 0, int ncomp = 1, int finest_level = -1) const;

};

#endif /*_NyxParticleContainer_H_*/
