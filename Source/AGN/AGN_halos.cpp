#ifdef AGN
#include <iomanip>
#include <algorithm>
#include <vector>
#include <iostream>
#include <string>
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
#include <Derive.H>

#include <AGNParticleContainer.H>

#if BL_USE_MPI
#include <MemInfo.H>
#endif

#ifdef REEBER
#ifdef REEBER_HIST
#include <ReeberAnalysis.H>
#endif /* REEBER */
#include <reeber.H>
#endif /* REEBER_HIST */

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
      if (comp != Density_comp)
        {
          MultiFab::Divide(state, state, Density_comp, comp, ncomp1, nghost);
        }
    }
}

void
Nyx::primitive_to_conserved(amrex::MultiFab& state)
{
   BL_PROFILE("Nyx::primitive_to_conserved()");
   int nghost = state.nGrow();

   // Multiply every component but Density by Density
   for (int comp = 0; comp < state.nComp(); comp++)
   {
      if (comp != Density_comp)
        {
          MultiFab::Multiply(state, state, Density_comp, comp, ncomp1, nghost);
        }
   }
}

void
Nyx::agn_halo_find (Real dt)
{
   BL_PROFILE("Nyx::agn_halo_find()");

   std::cerr << "ACHTUNG halo_find" << std::endl;

   const Real * dx = geom.CellSize();

   amrex::MultiFab& new_state = get_new_data(State_Type);
   const BoxArray& simBA = new_state.boxArray();
   const DistributionMapping& simDM = new_state.DistributionMap();
   int simComp = new_state.nComp();

   // First copy the existing state into orig_state.
   MultiFab orig_state(simBA, simDM, simComp, nghost1);
   MultiFab::Copy(orig_state, new_state,
                  comp0, comp0, simComp, nghost1);

   // These are passed into the AGN particles' Redistribute
   int lev_min = 0;
   int lev_max = 0;
   int ngrow   = 0;

#ifdef REEBER
#ifndef REEBER_HIST

   bool created_file = false;
   bool do_analysis = true;

   if (do_analysis || (reeber_int > 0 && nStep() % reeber_int == 0))
   {


       //       amrex::MultiFab reeberMF(simBA, simDM, 8, 0);
     //       amrex::MultiFab
     Vector<std::unique_ptr<MultiFab> > particle_mf;
     //     std::unique_ptr<MultiFab> particle_mf(new MultiFab(simBA, simDM, 4, 0));
       int cnt = 1;

       Real cur_time = state[State_Type].curTime();

       std::vector<Halo> reeber_halos;

       // global values updated when compute_average_density() called:
       /*
       average_gas_density;
       average_dm_density;
       average_neutr_density;
       average_total_density;
       */

       amrex::Vector<amrex::MultiFab*> state_levels;
       // Derive quantities and store in components 1... of MultiFAB

       amrex::Vector<amrex::IntVect> level_refinements;

       state_levels.resize(parent->finestLevel()+1);
       //Write all levels to Vector
       for (int lev = 0; lev <= parent->finestLevel(); lev++)
           {
             state_levels[lev]=&((get_level(lev)).get_new_data(State_Type));
         level_refinements.push_back(get_level(lev).fineRatio());
           }

       //This will fill all levels in particle_mf with lev_min=0
       Nyx::theDMPC()->AssignDensity(particle_mf, 0, 4, -1, 0);

       runReeberAnalysis(state_levels, particle_mf, Geom(), level_refinements, nStep(), do_analysis, reeber_halos);
#else

   const auto& reeber_density_var_list = getReeberHaloDensityVars();
   bool do_analysis(doAnalysisNow());

   bool created_file = false;

   if (do_analysis || (reeber_int > 0 && nStep() % reeber_int == 0))
   {

     // Before creating new AGN particles, check if any of the existing AGN particles should be merged
     agn_halo_merge();

     // Before creating new AGN particles, accrete mass onto existing particles
     agn_halo_accrete(dt);


       BoxArray reeberBA;
       DistributionMapping reeberDM;
       getAnalysisDecomposition(Geom(), ParallelDescriptor::NProcs(), reeberBA, reeberDM);
       cout << "getAnalysisDecomposition returns BoxArray size " << reeberBA.size() << endl;
       cout << "getAnalysisDecomposition returns DistributionMapping ProcessorMap size " << reeberDM.ProcessorMap().size() << endl;
       amrex::MultiFab reeberMF(reeberBA, reeberDM, reeber_density_var_list.size() + 1, nghost0);
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

#endif /* REEBER_HIST */
#else

       // Before creating new AGN particles, check if any of the existing AGN particles should be merged
       agn_halo_merge();

       // Before creating new AGN particles,
       // accrete mass and momentum onto existing particles.
       // No change to state.
       agn_halo_accrete(dt);

       // Here we just create place-holders for the halos which should come from REEBER
       std::vector<IntVect> reeber_halos_pos;
       std::vector<Real>    reeber_halos_mass;

       Real haloMass = 1.1e11;
       const Box& domainBox = grids.minimalBox();
       const IntVect& lo = domainBox.smallEnd();
       const IntVect& hi = domainBox.bigEnd();
       // With box 0:31, will have 0 and 31/2 = 15 and 31.
       int numSplits = 2;
       IntVect stepVect = (hi - lo) / numSplits;
       std::vector<IntVect> vertices(numSplits+1);
       for (int i = 0; i < numSplits; i++)
         {
           vertices[i] = lo + i * stepVect;
         }
       vertices[numSplits] = hi;
       Box vertBox(IntVect::Zero, numSplits*IntVect::Unit);
       for (BoxIterator bit(vertBox); bit.ok(); ++bit)
         {
           IntVect vert = bit();
           IntVect iv(D_DECL(vertices[vert[0]][0],
                             vertices[vert[1]][1],
                             vertices[vert[2]][2]));
           reeber_halos_pos.push_back(iv);
           reeber_halos_mass.push_back(haloMass);
         }

#endif // ifdef REEBER

       amrex::Real    halo_mass;
       amrex::IntVect halo_pos ;

       std::ofstream os;

       //std::cout << mass_halo_min << '\t' << mass_seed << std::endl;
       //std::cout << " *************************************** " << std::endl;

       // agn_density_old will hold the density from depositing the
       // mass of existing particles.
       MultiFab agn_density_old(simBA, simDM, ncomp1, nghost1);
       agn_density_old.setVal(0.0);

       // Deposit the mass now in the particles onto agn_density_old, on grid.
       // (No change to mass of particles.)
       Nyx::theAPC()->AssignDensitySingleLevel(agn_density_old, level);

       // Make sure the density put into ghost cells is added to valid regions
       agn_density_old.SumBoundary(geom.periodicity());

       // Convert new_state to primitive variables: rho, velocity, energy/rho.
       conserved_to_primitive(new_state);

       // Add agn_density_old to new_state, which holds primitive variables.
       // This is from depositing mass of existing particles.
       // Later, we'll subtract the deposited mass of all particles, old & new.
       amrex::MultiFab::Add(new_state, agn_density_old,
                            comp0, Density_comp, ncomp1, nghost0);

#ifdef REEBER
       for (const Halo& h : reeber_halos)
       {
#if 0
           // We aren't actually writing to this file so don't create it
           if (reeber_halos_pos.size() > 0)
           {
              if (!created_file)
                 os.open(amrex::Concatenate(amrex::Concatenate("debug-halos-", nStep(), 5), ParallelDescriptor::MyProc(), 2));
              created_file = true;
           }
#endif
           halo_mass = h.total_mass;
           halo_pos  = h.position;
#else


       for (int i = 0; i < reeber_halos_pos.size(); i++)
       {
           halo_mass = reeber_halos_mass[i];
           halo_pos  = reeber_halos_pos[i];
#endif

           if (halo_mass > mass_halo_min)
           {
                amrex::Real x = (halo_pos[0]+0.5) * dx[0];
                amrex::Real y = (halo_pos[1]+0.5) * dx[1];
                amrex::Real z = (halo_pos[2]+0.5) * dx[2];

                amrex::Real mass = mass_seed;

                int lev = 0;
                int grid = 0;
                int tile = 0;

                // Note that we are going to add the particle into grid 0 and tile 0 at level 0 --
                //      this is not actually where the particle belongs, but we will let the Redistribute call
                //      put it in the right place

                Nyx::theAPC()->AddOneParticle(lev,grid,tile,mass,x,y,z); // ,u,v,w);
                //std::cout << "ADDED A PARTICLE AT " << x << " " << y << " " << z << " WITH MASS " << mass << std::endl;
           }
       } // end of loop over creating new particles from halos

       //std::cout << " *************************************** " << std::endl;
       //std::cout << "  " << std::endl;

       // At this point the particles have all been created on the same process as the halo they came from,
       // but they are not on the "right" process for going forward

       // Call Redistribute so that the new particles get their cell, grid and process defined
       Nyx::theAPC()->Redistribute(lev_min, lev_max, ngrow);

       // Fill the "ghosts" vector with particles in ghost cells of each grid
       Nyx::theAPC()->fillNeighbors();

       // ComputeOverlap sets the ID of a particle to -1 if it is less than "cutoff" away from another
       //   particle and if it is newer than that particle
       Nyx::theAPC()->ComputeOverlap(level);

       // Clear the Neighbor Particle data structure
       Nyx::theAPC()->clearNeighbors();

       // This Redistribute is used to remove particles whose ID's have been set to -1 in ComputeOverlap
       Nyx::theAPC()->Redistribute(lev_min, lev_max, ngrow);

       // agn_density will hold the density we're going to remove from the grid.
       MultiFab agn_density(simBA, simDM, ncomp1, nghost1);
       agn_density.setVal(0.0);

       // Deposit the mass now in the particles onto the grid.
       // (No change to mass of particles.)
       Nyx::theAPC()->AssignDensitySingleLevel(agn_density, level);

       // Make sure the density put into ghost cells is added to valid regions
       agn_density.SumBoundary(geom.periodicity());

       // Take away the density from the gas that was added to the AGN particle:
       // density is in new_state, which holds primitive variables.
       amrex::MultiFab::Subtract(new_state, agn_density,
                                 comp0, Density_comp, ncomp1, nghost0);

       // Convert new_state to conserved variables: rho, momentum, energy.
       primitive_to_conserved(new_state);

       //Print() << "Going into ComputeParticleVelocity (no energy), number of AGN particles on this proc is "
       //       << Nyx::theAPC()->TotalNumberOfParticles(true, true) << endl;

       // Re-set the particle velocity (but not energy) after accretion,
       // using change of momentum density from orig_state to new_state,
       // which hold conserved variables.
       // No change to state, other than filling ghost cells.
       int add_energy = 0;
       Nyx::theAPC()->ComputeParticleVelocity(level, orig_state, new_state, add_energy);

       //Print() << "Going into ReleaseEnergy, number of AGN particles on this proc is "
       //       << Nyx::theAPC()->TotalNumberOfParticles(true, true) << endl;

       // AGN particles: may zero out energy.
       // new_state: may increase internal and total energy.
       MultiFab& D_new = get_new_data(DiagEOS_Type);
       Real a = get_comoving_a(new_a_time);
       Nyx::theAPC()->ReleaseEnergy(level, new_state, D_new, a);
       // Now new_state = get_new_data(State_Type) has been updated.

       //       Print() << "Step " << nStep() << " at end of Nyx_halos:" << endl;
       //       Nyx::theAPC()->writeAllAtLevel(level);

#ifdef REEBER

     }
#endif // ifdef REEBER
}

void
Nyx::agn_halo_merge ()
{
   int npart = Nyx::theAPC()->TotalNumberOfParticles(true, true);
   Nyx::theAPC()->fillNeighbors();
   Nyx::theAPC()->Merge(level);
   Nyx::theAPC()->clearNeighbors();

   // Call Redistribute to remove any particles with id = -1 (as set inside the Merge call)
   Nyx::theAPC()->Redistribute(level, level, nghost0);

   int nmerge = npart -Nyx::theAPC()->TotalNumberOfParticles(true, true);
   Print() << nStep() << '\t' << nmerge << endl;

}

void
Nyx::agn_halo_accrete (Real dt)
{
   amrex::MultiFab& new_state = get_new_data(State_Type);
   const BoxArray& simBA = new_state.boxArray();
   const DistributionMapping& simDM = new_state.DistributionMap();
   int ncomp = new_state.nComp();

   // First copy the existing state (new_state) into orig_state.
   MultiFab orig_state(simBA, simDM, ncomp, nghost1);
   MultiFab::Copy(orig_state, new_state,
                  comp0, comp0, ncomp, nghost1);

   // Convert new_state to primitive variables: rho, velocity, energy/rho.
   conserved_to_primitive(new_state);

   // Create a MultiFab to hold the density we're going to remove from the grid
   MultiFab agn_density_lost(simBA, simDM, ncomp1, nghost1);
   agn_density_lost.setVal(0.0);

   //Print() << "Going into AccreteMass, number of AGN particles on this proc is "
   //       << Nyx::theAPC()->TotalNumberOfParticles(true, true) << endl;

   // AGN particles: increase mass and energy.
   // new_state: no change, other than filling in ghost cells.
   // agn_density_lost: gets filled in.
   Nyx::theAPC()->AccreteMass(level, new_state, agn_density_lost, dt);

   // Make sure the density put into ghost cells is added to valid regions
   agn_density_lost.SumBoundary(geom.periodicity());

   // Take away the density from the gas that was added to the AGN particle.
   amrex::MultiFab::Subtract(new_state, agn_density_lost,
                             comp0, Density_comp, ncomp1, nghost0);

   // Convert new_state to conserved variables: rho, momentum, energy.
   primitive_to_conserved(new_state);

   // Re-set the particle velocity and energy after accretion,
   // using change of momentum density in state.
   // No change to state, other than filling ghost cells.
   int add_energy = 1;
   //Print() << "Going into ComputeParticleVelocity (and energy), number of AGN particles on this proc is "
   //       << Nyx::theAPC()->TotalNumberOfParticles(true, true) << endl;
   Nyx::theAPC()->ComputeParticleVelocity(level, orig_state, new_state, add_energy);
   // Now new_state = get_new_data(State_Type) has been updated.
}
#endif // AGN
