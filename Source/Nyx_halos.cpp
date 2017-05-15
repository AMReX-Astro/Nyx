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

using namespace amrex;

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
       amrex::MultiFab reeberMF(ba, dm, reeber_density_var_list.size() + 1, 0);
       int cnt = 1;

       Real cur_time = state[State_Type].curTime();

       // Derive quantities and store in components 1... of MultiFAB
       for (auto it = reeber_density_var_list.begin(); it != reeber_density_var_list.end(); ++it)
       {
           std::unique_ptr<MultiFab> derive_dat = particle_derive(*it, cur_time, 0);
           reeberMF.copy(*derive_dat, 0, cnt, 1, 0, 0);
           cnt++;
       }

       reeberMF.setVal(0, 0, 1, 0);
       for (int comp = 1; comp < reeberMF.nComp(); ++comp)
           amrex::MultiFab::Add(reeberMF, reeberMF, comp, 0, 1, 0);

       std::vector<Halo> reeber_halos;
       runReeberAnalysis(reeberMF, Geom(), nStep(), do_analysis, &reeber_halos);

#else

       // Before creating new AGN particles, check if any of the existing AGN particles should be merged
       halo_merge();

       cout << "Before accrete :" << endl;
       Nyx::theAPC()->writeAllAtLevel(level);

       // Before creating new AGN particles, accrete mass and momentum onto existing particles 
       halo_accrete(dt);

       cout << "After accrete :" << endl;
       Nyx::theAPC()->writeAllAtLevel(level);

       // Here we just create place-holders for the halos which should come from REEBER
       int num_halos = 10;
       std::vector<IntVect> reeber_halos_pos(num_halos);
       std::vector<Real>    reeber_halos_mass(num_halos);

       int i = 0;
       for (IntVect& iv : reeber_halos_pos)
       {
            i++;
            iv = IntVect(i+1,2*i+1,i+16);
       }

       for (Real& m : reeber_halos_mass)
            m = 1.1e11;

#endif // ifdef REEBER

       amrex::Real    halo_mass;
       amrex::IntVect halo_pos ;

       std::ofstream os;

       MultiFab new_state(simBA, simDM, simMF.nComp(), 1);
       MultiFab::Copy(new_state,simMF,0,0,simMF.nComp(),1);

       // Divide all components of new_state, other than Density, by density.
       for (int comp = 0; comp < new_state.nComp(); comp++)
         {
           if (comp != Density)
             {
               MultiFab::Divide(new_state, new_state, Density, comp, 1, 1);
             }
         }

       // Create a MultiFab to hold the density we're going to remove from the grid
       MultiFab agn_density(simBA, simDM, 1, 1);
       agn_density.setVal(0.0);

       // Deposit the mass now in the particles onto the grid (this doesn't change the mass of the particles)
       Nyx::theAPC()->AssignDensitySingleLevel(agn_density, 0);

       // Make sure the density put into ghost cells is added to valid regions
       agn_density.SumBoundary();

       // Add the density from the gas that is currently in the AGN particles.
       amrex::MultiFab::Add(new_state,agn_density,0,Density,1,0);

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
       Nyx::theAPC()->Redistribute(lev_min,lev_max,ngrow);

       Nyx::theAPC()->fillGhosts(level);
       Nyx::theAPC()->ComputeOverlap(level);
       Nyx::theAPC()->clearGhosts(level);

       Nyx::theAPC()->Redistribute(lev_min,lev_max,ngrow);

       // Zero this out again
       agn_density.setVal(0.0);

       // Deposit the mass now in the particles onto the grid (this doesn't change the mass of the particles)
       Nyx::theAPC()->AssignDensitySingleLevel(agn_density, 0);

       // Make sure the density put into ghost cells is added to valid regions
       agn_density.SumBoundary();

       // Take away the density from the gas that was added to the AGN particle.
       amrex::MultiFab::Subtract(new_state,agn_density,0,Density,1,0);

       // In new_state, everything but Density was divided by Density
       // at the beginning, and then Density was changed.
       // Now multiply everything but Density by Density.
       for (int comp = 0; comp < new_state.nComp(); comp++)
         {
           if (comp != Density)
             {
               MultiFab::Multiply(new_state, new_state, Density, comp, 1, 1);
             }
         }
       
       int add_energy = 0;
       Nyx::theAPC()->ComputeParticleVelocity(level,simMF,new_state,add_energy);

       MultiFab::Copy(simMF, new_state, 0, 0, simMF.nComp(), 0);

       cout << "At End of Nyx_halos:" << endl;
       Nyx::theAPC()->writeAllAtLevel(level);

       const amrex::Real time2 = ParallelDescriptor::second();
       if (ParallelDescriptor::IOProcessor())
         std::cout << std::endl << "===== Time to post-process: " << time2 - time1 << " sec" << std::endl;

#ifdef REEBER
     } else { // we have sidecars, so do everything in-transit

       int sidecarSignal(NyxHaloFinderSignal);
       const int MPI_IntraGroup_Broadcast_Rank = ParallelDescriptor::IOProcessor() ? MPI_ROOT : MPI_PROC_NULL;
       ParallelDescriptor::Bcast(&sidecarSignal, 1, MPI_IntraGroup_Broadcast_Rank,
                                 ParallelDescriptor::CommunicatorInter(whichSidecar));

       Geometry geom(Geom());
       Geometry::SendGeometryToSidecar(&geom, whichSidecar);

       amrex::MultiFab reeberMF(grids, reeber_density_var_list.size(), 0);
       int cnt = 0;
       // Derive quantities and store in components 1... of MultiFAB
       for (auto it = reeber_density_var_list.begin(); it != reeber_density_var_list.end(); ++it)
       {
           std::unique_ptr<MultiFab> derive_dat = particle_derive(*it, cur_time, 0);
           reeberMF.copy(*derive_dat, 0, cnt, 1, 0, 0);
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
   Nyx::theAPC()->Redistribute(0,0,0);
}

void
Nyx::halo_accrete (Real dt)
{
   amrex::MultiFab& orig_state = get_new_data(State_Type);
   const BoxArray& origBA = orig_state.boxArray();
   const DistributionMapping& origDM = orig_state.DistributionMap();

   // First copy the existing state into new_state
   MultiFab new_state(origBA,origDM,orig_state.nComp(),1);
   MultiFab::Copy(new_state,orig_state,0,0,orig_state.nComp(),1);

   // Divide all components of new_state, other than Density, by density.
   for (int comp = 0; comp < new_state.nComp(); comp++)
     {
       if (comp != Density)
         {
           MultiFab::Divide(new_state, new_state, Density, comp, 1, 1);
         }
     }

   // Create a MultiFab to hold the density we're going to remove from the grid
   MultiFab agn_density(origBA,origDM,1,1);
   agn_density.setVal(0.0);

   // Deposit the mass now in the particles onto the grid (this doesn't change the mass of the particles)
   Nyx::theAPC()->AssignDensitySingleLevel(agn_density, 0);

   // Make sure the density put into ghost cells is added to valid regions
   agn_density.SumBoundary();

   // Add the density from the gas that is currently in the AGN particles.
   amrex::MultiFab::Add(new_state,agn_density,0,Density,1,0);

   // Increase the mass of existing particles
   Real eps_rad = 0.1;
   Nyx::theAPC()->AccreteMass(level,orig_state,eps_rad,dt);

   // Zero this out again
   agn_density.setVal(0.0);

   // Deposit the mass now in the particles onto the grid (this doesn't change the mass of the particles)
   Nyx::theAPC()->AssignDensitySingleLevel(agn_density, 0);

   // Multiply this by 1/(1-eps) since we remove Mdot*dt from the gas but we only added Mdot*dt*(1-eps) to the particle 
   Real fac = 1. / (1. - eps_rad);
   agn_density.mult(fac,0,1,1);

   // Make sure the density put into ghost cells is added to valid regions
   agn_density.SumBoundary();

   // Take away the density from the gas that was added to the AGN particle (recall the (1-eps) weighting).
   amrex::MultiFab::Subtract(new_state,agn_density,0,Density,1,0);

   // In new_state, everything but Density was divided by Density
   // at the beginning, and then Density was changed.
   // Now multiply everything but Density by Density.
   for (int comp = 0; comp < new_state.nComp(); comp++)
     {
       if (comp != Density)
         {
           MultiFab::Multiply(new_state, new_state, Density, comp, 1, 1);
         }
     }

   // Re-set the particle velocity after accretion
   int add_energy = 1;
   Nyx::theAPC()->ComputeParticleVelocity(level,orig_state,new_state, add_energy);
}
#endif // AGN
