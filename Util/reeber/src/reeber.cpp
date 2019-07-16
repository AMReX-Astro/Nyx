#include "Nyx.H"
#include "reeber.H"

// Reeber and DIY includes
#include "reeber-real.h"


#include <diy/master.hpp>
#include <diy/io/block.hpp>
#include <diy/io/shared.hpp>
#include <diy/decomposition.hpp>

#include <dlog/stats.h>
#include <dlog/log.h>
#include <opts/opts.h>
#include <error.h>
#include <AMReX_Geometry.H>

#include "fab-block.h"
#include "fab-cc-block.h"
#include "reader-interfaces.h"

#include "amr-connected-components-complex.h"


// block-independent types
using AMRLink = diy::AMRLink;

using Bounds = diy::DiscreteBounds;
using AmrVertexId = r::AmrVertexId;
using AmrEdge = reeber::AmrEdge;

#define DIM 3

using FabBlockR = FabBlock<Real, DIM>;

using Block = FabComponentBlock<Real, DIM>;
using Vertex = Block::Vertex;
using Component = Block::Component;
using MaskedBox = Block::MaskedBox;
using GidVector = Block::GidVector;
using GidSet = Block::GidSet;

using TripletMergeTree = Block::TripletMergeTree;
using Neighbor = TripletMergeTree::Neighbor;

struct IsAmrVertexLocal
{
    bool operator()(const Block& b, const Neighbor& from) const
    {
        return from->vertex.gid == b.gid;
    }
};

void compute_halos(diy::mpi::communicator& world,
        diy::Master& master_reader,
        int threads,
        diy::DiscreteBounds domain,
        Real absolute_rho,
        bool negate,
        Real min_halo_volume)
{
    world.barrier();
    std::string prefix = "./DIY.XXXXXX";
    int in_memory = -1;
    std::string log_level = "info";
    diy::FileStorage storage(prefix);

    dlog::add_stream(std::cerr, dlog::severity(log_level))
            << dlog::stamp() << dlog::aux_reporter(world.rank()) << dlog::color_pre() << dlog::level()
            << dlog::color_post() >> dlog::flush();


    // copy FabBlocks to FabComponentBlocks
    // in FabTmtConstructor mask will be set and local trees will be computed
    // FabBlock can be safely discarded afterwards
    diy::Master master(world, threads, in_memory, &Block::create, &Block::destroy, &storage, &Block::save,
            &Block::load);
    master_reader.foreach(
            [&master, domain, absolute_rho, negate](FabBlockR* b, const diy::Master::ProxyWithLink& cp) {
                auto* l = static_cast<AMRLink*>(cp.link());
                AMRLink* new_link = new AMRLink(*l);

                // prepare neighbor box info to save in MaskedBox
                // TODO: refinment vector
                int local_ref = l->refinement()[0];
                int local_lev = l->level();

                master.add(cp.gid(),
                        new Block(b->fab, b->extra_names_, b->extra_fabs_, local_ref, local_lev, domain,
                                l->bounds(),
                                l->core(), cp.gid(),
                                new_link, absolute_rho, negate, /*absolute = */ true),
                        new_link);

            });

    int global_n_undone = 1;

    master.foreach(&send_edges_to_neighbors_cc<Real, DIM>);
    master.exchange();
    master.foreach(&delete_low_edges_cc<Real, DIM>);

    int rounds = 0;
    while(global_n_undone)
    {
        rounds++;

        master.foreach(&amr_cc_send<Real, DIM>);
        master.exchange();
        master.foreach(&amr_cc_receive<Real, DIM>);
        master.exchange();

        global_n_undone = master.proxy(master.loaded_block()).read<int>();
    }

    master.foreach([](Block* b, const diy::Master::ProxyWithLink& cp) {
        b->compute_final_connected_components();
        b->compute_local_integral();
    });

    LOG_SEV_IF(world.rank() == 0, info) << "Local integrals computed";
    dlog::flush();

    bool has_density = true;
    bool has_particle_mass_density = true;
    bool has_momentum = false;

    master.foreach(
            [domain, min_halo_volume,
                    has_density, has_particle_mass_density,
                    has_momentum](
                    Block* b,
                    const diy::Master::ProxyWithLink& cp) {

                diy::Point<int, 3> domain_shape;
                for(int i = 0; i < 3; ++i)
                {
                    domain_shape[i] = domain.max[i] - domain.min[i] + 1;
                }

                diy::GridRef<void*, 3> domain_box(nullptr, domain_shape, /* c_order = */ false);

                // local integral already stores number of vertices (set in init)
                // so we add it here just to print it
                b->extra_names_.insert(b->extra_names_.begin(), std::string("n_vertices"));

                for(const auto& root_values_pair : b->local_integral_)
                {
                    AmrVertexId root = root_values_pair.first;
                    if (root.gid != b->gid)
                        continue;

                    auto& values = root_values_pair.second;

                    Real n_vertices = values.at("n_vertices");
                    Real halo_volume = values.at("n_vertices_sf");

                    if (halo_volume < min_halo_volume)
                        continue;

                    Real m_gas = has_density ? values.at("density") : 0;
                    Real m_particles = has_particle_mass_density ? values.at("particle_mass_density") : 0;
                    Real m_total = m_gas + m_particles;
                    // TODO: add halo to vector

                }
            });
}


using namespace amrex;

/* // Global components of state
int Nyx::Density = -1;
int Nyx::Eden = -1;
int Nyx::Eint = -1;
int Nyx::Xmom = -1;
int Nyx::Ymom = -1;
int Nyx::Zmom = -1;*/

/* // Global density values
   average_gas_density;
   average_dm_density;
   average_neutr_density;
   average_total_density;*/

void Nyx::runReeberAnalysis(Vector<MultiFab*>& new_state, Vector<std::unique_ptr<MultiFab> >& particle_mf,
        const amrex::Geometry geom, int n_step, bool do_analysis, std::vector<Halo>& reeber_halos)
{

  if(verbose)
    amrex::Print()<<"Running Reeber anaylsis"<<std::endl;
  for (int lev = 0; lev <= parent->finestLevel(); lev++)
    {
      for ( MFIter mfi(*(new_state[lev]), TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {

      const Box& tbx = mfi.tilebox();
      Array4<const Real> state = new_state[lev]->array(mfi);
      Array4<Real>  particle = particle_mf[lev]->array(mfi);
      /*
          FArrayBox& fab_state = new_state[lev][mfi];
          FArrayBox& fab_particle = particle_mf[lev][mfi];*/
      const Dim3 lo = amrex::lbound(tbx);
      const Dim3 hi = amrex::ubound(tbx);
      int ncomp = 4;

      for (int n = 0; n < ncomp; ++n) {
        for (int z = lo.z; z <= hi.z; ++z) {
          for (int y = lo.y; y <= hi.y; ++y) {
        AMREX_PRAGMA_SIMD
          for (int x = lo.x; x <= hi.x; ++x) {
            state(x,y,z,n);
            particle(x,y,z,n);
          }
          }
        }
      }
    }
    }
  return;
}
