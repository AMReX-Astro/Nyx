#ifdef AGN
#include <iomanip>
#include <algorithm>
#include <vector>
#include <iostream>
#include <string>
#include <unistd.h>
#include <math.h>

using std::cout;
using std::cerr;
using std::endl;
using std::istream;
using std::ostream;
using std::pair;
using std::string;

#include <AMReX_CONSTANTS.H>
#include <Nyx.H>
#include <Nyx_F.H>
#include <Derive_F.H>

#include "AGNParticleContainer.H"

#if BL_USE_MPI
#include "MemInfo.H"
#endif

#ifdef REEBER
#include <ReeberAnalysis.H>
#endif /* REEBER */

const int NyxHaloFinderSignal = 42;

// For debugging.
const int nghost0 = 0;
const int nghost1 = 1;
const int ncomp1 = 1;
const int comp0 = 0;

using namespace amrex;

void
Nyx::conserved_to_primitive(amrex::MultiFab& state)
{ // divide every component but Density by Density
  int nghost = state.nGrow();
  for (int comp = 0; comp < state.nComp(); comp++)
    {
      if (comp != Density)
        {
          MultiFab::Divide(state, state, Density, comp, ncomp1, nghost);
        }
    }
}

void
Nyx::primitive_to_conserved(amrex::MultiFab& state)
{ // multiply every component but Density by Density
  int nghost = state.nGrow();
  for (int comp = 0; comp < state.nComp(); comp++)
    {
      if (comp != Density)
        {
          MultiFab::Multiply(state, state, Density, comp, ncomp1, nghost);
        }
    }
}

void
Nyx::halo_find (Real dt)
{
   BL_PROFILE("Nyx::halo_find()");

   const int whichSidecar(0);

   const amrex::Real time1 = ParallelDescriptor::second();

   const Real * dx = geom.CellSize();

   amrex::MultiFab& simMF = get_new_data(State_Type);
   const BoxArray& simBA = simMF.boxArray();
   const DistributionMapping& simDM = simMF.DistributionMap();

   // These are passed into the AGN particles' Redistribute
   int lev_min = 0;
   int lev_max = 0;
   int ngrow   = 0;

#ifdef REEBER
   const auto& reeber_density_var_list = getReeberHaloDensityVars();
   bool do_analysis(doAnalysisNow());

   bool created_file = false;

   if (do_analysis || (reeber_int > 0 && nStep() % reeber_int == 0)) 
   {

     // Before creating new AGN particles, check if any of the existing AGN particles should be merged
     halo_merge();

     // Before creating new AGN particles, accrete mass onto existing particles 
     halo_accrete(dt);

     if (ParallelDescriptor::NProcsSidecar(0) <= 0) 
     { // we have no sidecars, so do everything in situ

       BoxArray ba;
       DistributionMapping dm;
       getAnalysisDecomposition(Geom(), ParallelDescriptor::NProcs(), ba, dm);
       cout << "getAnalysisDecomposition returns BoxArray size " << ba.size() << endl;
       cout << "getAnalysisDecomposition returns DistributionMapping ProcessorMap size " << dm.ProcessorMap().size() << endl;
       amrex::MultiFab reeberMF(ba, dm, reeber_density_var_list.size() + 1, nghost0);
       int cnt = 1;

       Real cur_time = state[State_Type].curTime();

       // Derive quantities and store in components 1... of MultiFAB
       for (auto it = reeber_density_var_list.begin(); it != reeber_density_var_list.end(); ++it)
       {
           std::unique_ptr<MultiFab> derive_dat = particle_derive(*it, cur_time, 0);
           reeberMF.copy(*derive_dat, comp0, cnt, ncomp1, nghost0, nghost0);
           cnt++;
       }

       reeberMF.setVal(0.0, comp0, ncomp1, nghost0);
       for (int comp = 1; comp < reeberMF.nComp(); ++comp)
         {
           amrex::MultiFab::Add(reeberMF, reeberMF,
                                comp, comp0, ncomp1, nghost0);
         }

       std::vector<Halo> reeber_halos;
       runReeberAnalysis(reeberMF, Geom(), nStep(), do_analysis, &reeber_halos);

#else

       // Before creating new AGN particles, check if any of the existing AGN particles should be merged
       halo_merge();

       cout << "Before accrete :" << endl;
       Nyx::theAPC()->writeAllAtLevel(level);

       // Before creating new AGN particles,
       // accrete mass and momentum onto existing particles.
       // No change to state.
       halo_accrete(dt);

       cout << "After accrete :" << endl;
       Nyx::theAPC()->writeAllAtLevel(level);

       // Here we just create place-holders for the halos which should come from REEBER
       std::vector<IntVect> reeber_halos_pos;
       std::vector<Real>    reeber_halos_mass;

       Real haloMass = 1.1e11;
       const Box& domainBox = grids.minimalBox();
       const IntVect& lo = domainBox.smallEnd();
       const IntVect& hi = domainBox.bigEnd();
       Box cornerBox(IntVect::Zero, IntVect::Unit);
       for (BoxIterator bit(cornerBox); bit.ok(); ++bit)
         {
           IntVect iv = lo + (hi - lo) * bit();
           reeber_halos_pos.push_back(iv);
           reeber_halos_mass.push_back(haloMass);
         }

#endif // ifdef REEBER

       amrex::Real    halo_mass;
       amrex::IntVect halo_pos ;

       std::ofstream os;

       MultiFab new_state(simBA, simDM, simMF.nComp(), nghost1);
       MultiFab::Copy(new_state, simMF,
                      comp0, comp0, simMF.nComp(), nghost1);
       // Convert new_state to primitive variables: rho, velocity, energy/rho.
       conserved_to_primitive(new_state);

       std::cout << "  " << std::endl;
       std::cout << " *************************************** " << std::endl;

#ifdef REEBER
       for (const Halo& h : reeber_halos)
       {
           if (!created_file)
              os.open(amrex::Concatenate(amrex::Concatenate("debug-halos-", nStep(), 5), ParallelDescriptor::MyProc(), 2));
           created_file = true;
           halo_mass = h.totalMass;
           halo_pos  = h.position;
#else

       // Now loop over the halos
       for (int i = 0; i < reeber_halos_pos.size(); i++)
       {
           halo_mass = reeber_halos_mass[i];
           halo_pos  = reeber_halos_pos[i];
#endif

           if (halo_mass > 1.e10)
           {
                amrex::Real x = (halo_pos[0]+0.5) * dx[0];
                amrex::Real y = (halo_pos[1]+0.5) * dx[1];
                amrex::Real z = (halo_pos[2]+0.5) * dx[2];
   
                amrex::Real mass = 1.e5;

                int lev = 0;
                int grid = 0;
                int tile = 0;

                // Note that we are going to add the particle into grid 0 and tile 0 at level 0 -- 
                //      this is not actually where the particle belongs, but we will let the Redistribute call
                //      put it in the right place

                Nyx::theAPC()->AddOneParticle(lev,grid,tile,mass,x,y,z); // ,u,v,w);
                std::cout << "ADDED A PARTICLE AT " << x << " " << y << " " << z << " WITH MASS " << mass << std::endl;
           }
       } // end of loop over creating new particles from halos

       std::cout << " *************************************** " << std::endl;
       std::cout << "  " << std::endl;

       // At this point the particles have all been created on the same process as the halo they came from,
       // but they are not on the "right" process for going forward

       // Call Redistribute so that the new particles get their cell, grid and process defined
       Nyx::theAPC()->Redistribute(lev_min, lev_max, ngrow);

       Nyx::theAPC()->fillGhosts(level);
       Nyx::theAPC()->ComputeOverlap(level);
       Nyx::theAPC()->clearGhosts(level);

       Nyx::theAPC()->Redistribute(lev_min, lev_max, ngrow);

       // agn_density will hold the density we're going to remove from the grid.
       MultiFab agn_density(simBA, simDM, ncomp1, nghost1);
       agn_density.setVal(0.0);

       // Deposit the mass now in the particles onto the grid.
       // (No change to mass of particles.)
       Nyx::theAPC()->AssignDensitySingleLevel(agn_density, level);

       // Make sure the density put into ghost cells is added to valid regions
       agn_density.SumBoundary(geom.periodicity());

       // Take away the density from the gas that was added to the AGN particle.
       amrex::MultiFab::Subtract(new_state, agn_density,
                                 comp0, Density, ncomp1, nghost0);

       // Convert new_state to conserved variables: rho, momentum, energy.
       // Since the density has changed, the other variables change accordingly.
       primitive_to_conserved(new_state);

       cout << "Going into ComputeParticleVelocity (no energy), number of AGN particles on this proc is "
            << Nyx::theAPC()->TotalNumberOfParticles(true, true) << endl;
       int add_energy = 0;
       Nyx::theAPC()->ComputeParticleVelocity(level, simMF, new_state, add_energy);

       Real T_min = 1.0e+7;
       cout << "Going into ReleaseEnergy, number of AGN particles on this proc is "
            << Nyx::theAPC()->TotalNumberOfParticles(true, true) << endl;
       // AGN particles: may zero out energy.
       // new_state: may increase internal and total energy.
       MultiFab& D_new = get_new_data(DiagEOS_Type);
       Real a = get_comoving_a(new_a_time);
       Nyx::theAPC()->ReleaseEnergy(level, new_state, D_new, a, T_min);

       MultiFab::Copy(simMF, new_state,
                      comp0, comp0, simMF.nComp(), nghost0);

       cout << "At End of Nyx_halos:" << endl;
       Nyx::theAPC()->writeAllAtLevel(level);

       const amrex::Real time2 = ParallelDescriptor::second();
       if (ParallelDescriptor::IOProcessor())
         std::cout << std::endl << "===== Time to post-process: " << time2 - time1 << " sec" << std::endl;

#ifdef REEBER
     } 
#if 0
     else { // we have sidecars, so do everything in-transit

       int sidecarSignal(NyxHaloFinderSignal);
       const int MPI_IntraGroup_Broadcast_Rank = ParallelDescriptor::IOProcessor() ? MPI_ROOT : MPI_PROC_NULL;
       ParallelDescriptor::Bcast(&sidecarSignal, 1, MPI_IntraGroup_Broadcast_Rank,
                                 ParallelDescriptor::CommunicatorInter(whichSidecar));

       Geometry geom(Geom());
       Geometry::SendGeometryToSidecar(&geom, whichSidecar);

       // FIXME: What is distribution mapping?
       amrex::MultiFab reeberMF(grids, reeber_density_var_list.size(), 0);
       int cnt = 0;
       // Derive quantities and store in components 1... of MultiFAB
       for (auto it = reeber_density_var_list.begin(); it != reeber_density_var_list.end(); ++it)
       {
           std::unique_ptr<MultiFab> derive_dat = particle_derive(*it, cur_time, 0);
           reeberMF.copy(*derive_dat, comp0, cnt, ncomp1, nghost0, nghost0);
           cnt++;
       }

       int time_step(nStep()), nComp(reeberMF.nComp());

       ParallelDescriptor::Bcast(&nComp, 1, MPI_IntraGroup_Broadcast_Rank,
                                 ParallelDescriptor::CommunicatorInter(whichSidecar));

       amrex::MultiFab *mfSource = &reeberMF;
       amrex::MultiFab *mfDest = 0;
       int srcComp(0), destComp(1);
       int srcNGhost(0), destNGhost(0);
       MPI_Comm commSrc(ParallelDescriptor::CommunicatorComp());
       MPI_Comm commDest(ParallelDescriptor::CommunicatorSidecar());
       MPI_Comm commInter(ParallelDescriptor::CommunicatorInter(whichSidecar));
       MPI_Comm commBoth(ParallelDescriptor::CommunicatorBoth(whichSidecar));
       bool isSrc(true);

       amrex::MultiFab::copyInter(mfSource, mfDest, srcComp, destComp, nComp,
                           srcNGhost, destNGhost,
                           commSrc, commDest, commInter, commBoth,
                           isSrc);


       ParallelDescriptor::Bcast(&time_step, 1, MPI_IntraGroup_Broadcast_Rank,
                                 ParallelDescriptor::CommunicatorInter(whichSidecar));

       int do_analysis_bcast(do_analysis);
       ParallelDescriptor::Bcast(&do_analysis_bcast, 1, MPI_IntraGroup_Broadcast_Rank,
                                 ParallelDescriptor::CommunicatorInter(whichSidecar));

       amrex::Real time2(ParallelDescriptor::second());
       if(ParallelDescriptor::IOProcessor()) 
         std::cout << "COMPUTE PROCESSES: time spent sending data to sidecars: " << time2 - time1 << std::endl;

#endif // if 0
     }
#endif // ifdef REEBER
}

void
Nyx::halo_merge ()
{
   Nyx::theAPC()->fillGhosts(level);
   Nyx::theAPC()->Merge(level);
   Nyx::theAPC()->clearGhosts(level);

   // Call Redistribute to remove any particles with id = -1 (as set inside the Merge call)
   Nyx::theAPC()->Redistribute(level, level, nghost0);
}

void
Nyx::halo_accrete (Real dt)
{
   amrex::MultiFab& orig_state = get_new_data(State_Type);
   const BoxArray& origBA = orig_state.boxArray();
   const DistributionMapping& origDM = orig_state.DistributionMap();

   // First copy the existing state into new_state
   MultiFab new_state(origBA, origDM, orig_state.nComp(), nghost1);
   MultiFab::Copy(new_state, orig_state,
                  comp0, comp0, orig_state.nComp(), nghost1);

   // Convert new_state to primitive variables: rho, velocity, energy/rho.
   conserved_to_primitive(new_state);

   // Create a MultiFab to hold the density we're going to remove from the grid
   MultiFab agn_density_lost(origBA, origDM, ncomp1, nghost1);
   agn_density_lost.setVal(0.0);

   Real eps_rad = 0.1;
   Real eps_coupling = 0.15;
   cout << "Going into AccreteMass, number of AGN particles on this proc is "
        << Nyx::theAPC()->TotalNumberOfParticles(true, true) << endl;
   // AGN particles: increase mass and energy.
   // new_state: no change, other than filling in ghost cells.
   // agn_density_lost: gets filled in.
   Nyx::theAPC()->AccreteMass(level, new_state, agn_density_lost,
                              eps_rad, eps_coupling, dt);

   // Make sure the density put into ghost cells is added to valid regions
   agn_density_lost.SumBoundary(geom.periodicity());

   // Take away the density from the gas that was added to the AGN particle.
   amrex::MultiFab::Subtract(new_state, agn_density_lost,
                             comp0, Density, ncomp1, nghost0);

   // Convert new_state to conserved variables: rho, momentum, energy.
   primitive_to_conserved(new_state);

   // Re-set the particle velocity after accretion
   int add_energy = 1;
   cout << "Going into ComputeParticleVelocity (and energy), number of AGN particles on this proc is "
        << Nyx::theAPC()->TotalNumberOfParticles(true, true) << endl;
   Nyx::theAPC()->ComputeParticleVelocity(level, orig_state, new_state, add_energy);
}
#endif // AGN
