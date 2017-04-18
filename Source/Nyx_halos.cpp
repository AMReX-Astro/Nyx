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

#include "NyxParticleContainer.H"

#if BL_USE_MPI
#include "MemInfo.H"
#endif

#ifdef REEBER
#include <ReeberAnalysis.H>
#endif /* REEBER */

const int NyxHaloFinderSignal = 42;

using namespace amrex;

void
Nyx::halo_find ()
{
   BL_PROFILE("Nyx::halo_find()");

   const amrex::Real cur_time = state[State_Type].curTime();

   const int whichSidecar(0);

   const amrex::Real time1 = ParallelDescriptor::second();

   const Real * dx = geom.CellSize();

   const amrex::MultiFab& simMF = get_new_data(State_Type);
   const BoxArray& simBA = simMF.boxArray();
   const DistributionMapping& simDM = simMF.DistributionMap();

#ifdef REEBER
   const auto& reeber_density_var_list = getReeberHaloDensityVars();
   bool do_analysis(doAnalysisNow());

   if (do_analysis || (reeber_int > 0 && nStep() % reeber_int == 0)) {
     if (ParallelDescriptor::NProcsSidecar(0) <= 0) { // we have no sidecars, so do everything in situ

       BoxArray ba;
       DistributionMapping dm;
       getAnalysisDecomposition(Geom(), ParallelDescriptor::NProcs(), ba, dm);
       amrex::MultiFab reeberMF(ba, reeber_density_var_list.size() + 1, 0, dm);
       int cnt = 1;
       // Derive quantities and store in components 1... of MultiFAB
       for (auto it = reeber_density_var_list.begin(); it != reeber_density_var_list.end(); ++it)
       {
           amrex::MultiFab *derive_dat = particle_derive(*it, cur_time, 0); // FIXME: Is this the right way? 
           reeberMF.copy(*derive_dat, 0, cnt, 1, 0, 0);
           delete derive_dat;
           cnt++;
       }

       reeberMF.setVal(0, 0, 1, 0);
       for (int comp = 1; comp < reeberMF.nComp(); ++comp)
           amrex::MultiFab::Add(reeberMF, reeberMF, comp, 0, 1, 0);

       std::vector<Halo> reeber_halos;
       runReeberAnalysis(reeberMF, Geom(), nStep(), do_analysis, &reeber_halos);

       // Redistribute halos to "correct" processor for simulation
       // FIXME: This is a hack that maps things to amrex::MultiFabs and back. This should be
       // changed to use Nyx's particle redistribution infrastructure.
       reeberMF.setVal(0, 0, reeber_density_var_list.size() + 1, 0);
       // Deposit halo into amrex::MultiFab. This works because (i) There is only one box
       // per processors and halos returned by Reeber will always be in that box;
       // (ii) Halos have different positions; (iii) Halos have a mass that differs from
       // zero.
       BL_ASSERT(dm[ParallelDescriptor::MyProc()] == ParallelDescriptor::MyProc());
       FArrayBox& my_fab = reeberMF[ParallelDescriptor::MyProc()];

       for (const Halo& h : reeber_halos)
       {
           BL_ASSERT(reeberMF.fabbox(ParallelDescriptor::MyProc()).contains(h.position));
           my_fab(h.position, 0) = h.totalMass;
           for (int comp = 0; comp < reeber_density_var_list.size(); ++comp)
           {
               my_fab(h.position, comp + 1) = h.individualMasses[comp];
           }
       }

       // Actual redistribution
       //amrex::MultiFab redistributeFab(m_leveldata.boxArray(), reeber_density_var_list.size() + 1, 0, m_leveldata.DistributionMap());

       amrex::MultiFab redistributeFab(simBA, reeber_density_var_list.size() + 1, 0, simDM);
       redistributeFab.copy(reeberMF);
       // Re-extract halos
       reeber_halos.clear();
       for (MFIter mfi(redistributeFab); mfi.isValid(); ++mfi)
       {
           const Box& currBox = mfi.fabbox();
           for (IntVect iv = currBox.smallEnd(); iv <= currBox.bigEnd(); currBox.next(iv))
           {
               amrex::Real totalMass = redistributeFab[mfi](iv, 0);
               if (totalMass > 0)
               {
                   std::vector<amrex::Real> masses(reeber_density_var_list.size(), 0);
                   for (int comp = 0; comp < reeber_density_var_list.size(); ++comp)
                   { masses[comp] = redistributeFab[mfi](iv, comp + 1);
                   }
                   reeber_halos.emplace_back(iv, totalMass, masses);
               }
           }
        }
       // NOTE: ZARIJA, GET YOUR FRESH HALOS HERE!!!

#endif // ifdef REEBER

       amrex::Real    halo_mass;
       amrex::IntVect halo_pos ;

       amrex::Real mass, x, y, z;

       std::ofstream os;

#ifdef REEBER
       for (const Halo& h : reeber_halos)
       {
           if (!created_file)
              os.open(BoxLib::Concatenate(BoxLib::Concatenate("debug-halos-", nStep(), 5), ParallelDescriptor::MyProc(), 2));
           created_file = true;
           halo_mass = h.totalMass;
           halo_pos  = h.position;
#else
       {
           halo_mass = 1.1e11;
           halo_pos = IntVect(3,3,3);
#endif
           std::cout << "HALO HERE !!! " << halo_pos << " " << halo_mass << std::endl;

           bool new_particle_created = false; 

           if (halo_mass > 1.e10)
           {
                x = (halo_pos[0]+0.5) * dx[0];
                y = (halo_pos[1]+0.5) * dx[1];
                z = (halo_pos[2]+0.5) * dx[2];
   
                amrex::Real scaled_halo_mass = halo_mass/1.e13;
    
                mass = std::pow(10.0,8.18) * pow(scaled_halo_mass,1.55);

                int lev = 0;
                int grid = 0;
                int tile = 0;

                // Note that we are going to add the particle into grid 0 and tile 0 at level 0 -- 
                //      this is not actually where the particle belongs, but we will let the Redistribute call
                //      put it in the right place
                Nyx::theAPC()->AddOneParticle(lev,grid,tile,halo_mass,x,y,z); // ,u,v,w);

                std::cout << "  " << std::endl;
                std::cout << " *************************************** " << std::endl;
                std::cout << "ADDED A PARTICLE AT " << x << " " << y << " " << z << " WITH MASS " << mass << std::endl;
                std::cout << " *************************************** " << std::endl;
                std::cout << "  " << std::endl;

                new_particle_created = true; 
#if 0
              for (MFIter mfi(agn_density); mfi.isValid(); ++mfi)
              {
                  const Box& currBox = mfi.fabbox();
                  if (currBox.contains(halo_pos)) 
                  {
                     int index = mfi.index();

                     std::cout << " CELL / INDEX " << halo_pos << " " << index << std::endl;

                     // Figure out if there are already any nearby AGN particles 
                     amrex::Real agn_mass = 0.0;
                     for (int k  = halo_pos[2] - 1; k <= halo_pos[2] + 1; k++)
                      for (int j  = halo_pos[1] - 1; j <= halo_pos[1] + 1; j++)
                       for (int i  = halo_pos[0] - 1; i <= halo_pos[0] + 1; i++)
                       {
                         IntVect iv(i,j,k);
                         std::cout << "WHICH CELL " << iv << " " << agn_density[index](iv,0) << std::endl;  
                         agn_mass += agn_density[index](iv,0);
                       }

                     std::cout << "NEARBY AGN_MASS " << agn_mass   << std::endl;

                     if (!(agn_mass > 0.0))
                     {
                        gas_sum[0] = 0.0;
                        gas_sum[1] = 0.0;
                        gas_sum[2] = 0.0;
                        gas_sum[3] = 0.0;
                        for (int k  = halo_pos[2] - 1; k <= halo_pos[2] + 1; k++)
                         for (int j  = halo_pos[1] - 1; j <= halo_pos[1] + 1; j++)
                          for (int i  = halo_pos[0] - 1; i <= halo_pos[0] + 1; i++)
                          {
                            IntVect iv(i,j,k);
                            gas_sum[0] += gas_state[index](iv,Density);
                            gas_sum[1] += gas_state[index](iv,Xmom);
                            gas_sum[2] += gas_state[index](iv,Ymom);
                            gas_sum[3] += gas_state[index](iv,Zmom);
                          }
                        gas_sum[0] *= cellVol;
                        gas_sum[1] *= cellVol;
                        gas_sum[2] *= cellVol;
                        gas_sum[3] *= cellVol;

                        amrex::Real frac = mass / gas_sum[0];
                        std::cout << "FRAC OF GAS MASS " << frac << std::endl;

                        // Define the velocity of the new particle based on the momentum of the surrounding gas
                        u = gas_sum[1] / mass;
                        v = gas_sum[2] / mass;
                        w = gas_sum[3] / mass;

                        for (int k  = halo_pos[2] - 1; k <= halo_pos[2] + 1; k++)
                         for (int j  = halo_pos[1] - 1; j <= halo_pos[1] + 1; j++)
                          for (int i  = halo_pos[0] - 1; i <= halo_pos[0] + 1; i++)
                          {
                            IntVect iv(i,j,k);
                            sub_state[index](iv,0) -= frac * gas_state[index](iv,Density);
                            sub_state[index](iv,1) -= frac * gas_state[index](iv,Xmom);
                            sub_state[index](iv,2) -= frac * gas_state[index](iv,Ymom);
                            sub_state[index](iv,3) -= frac * gas_state[index](iv,Zmom);
                          }
                     }
                  } // End of test on whether particle is in this box
              } // End of MFIter
#endif
           } // End of test on halo mass
       } // End of loop over halos

       // Call Redistribute so that the new particles get their cell, grid and process defined
       Nyx::theAPC()->Redistribute(false,true,0);

#if 0
       // Take away the density/momentum from the gas that was added to the AGN particle.
       sub_state.SumBoundary();
       amrex::MultiFab::Add(gas_state,sub_state,0,Density,1+BL_SPACEDIM,0);

       // This is just a test to make sure the particle deposits mass
       MultiFab agn_density(simBA, simDM, 1, 1);
       agn_density.setVal(0.0);
       Nyx::theAPC()->AssignDensitySingleLevel(agn_density, 0);
       amrex::Real agn_mass_new = agn_density.norm0();

       if (ParallelDescriptor::IOProcessor())
          std::cout << "MAX NORM OF AGN_DENSITY AFTER HALO STUFF " << agn_mass_new << std::endl;
#endif

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
           amrex::MultiFab *derive_dat = particle_derive(*it, cur_time, 0); // FIXME: Is this the right way?
           reeberMF.copy(*derive_dat, 0, cnt, 1, 0, 0);
           delete derive_dat;
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
#endif // AGN
