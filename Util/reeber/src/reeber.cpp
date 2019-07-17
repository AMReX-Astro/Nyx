#include <vector>
#include <string>

#include "Nyx.H"
#include "reeber.H"

// Reeber and DIY includes
#include <reeber-real.h>


#include <diy/master.hpp>
#include <diy/io/block.hpp>
#include <diy/io/shared.hpp>
#include <diy/resolve.hpp>
#include <diy/decomposition.hpp>

#include <dlog/stats.h>
#include <dlog/log.h>
#include <opts/opts.h>
#include <error.h>
#include <AMReX_Geometry.H>

#include <fab-block.h>
#include <fab-cc-block.h>
//#include <reader-interfaces.h"

#include <amr-connected-components-complex.h>


using namespace amrex;

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


diy::AMRLink::Bounds bounds(const amrex::Box& box)
{
    diy::AMRLink::Bounds bounds(3);
    for(int i = 0; i < 3; ++i)
    {
        bounds.min[i] = box.loVect()[i];
        bounds.max[i] = box.hiVect()[i];
    }
    return bounds;
}

void record_wrap(const Box& domain, const Box& valid_box, const std::array<bool, DIM>& is_periodic, diy::AMRLink* link)
{
    for(int dir_x : {-1, 0, 1})
    {
        if (!is_periodic[0] && dir_x) continue;
        if (dir_x < 0 && valid_box.loVect()[0] != domain.loVect()[0]) continue;
        if (dir_x > 0 && valid_box.hiVect()[0] != domain.hiVect()[0]) continue;

        for(int dir_y : {-1, 0, 1})
        {
            if (!is_periodic[1] && dir_y) continue;
            if (dir_y < 0 && valid_box.loVect()[1] != domain.loVect()[1]) continue;
            if (dir_y > 0 && valid_box.hiVect()[1] != domain.hiVect()[1]) continue;
            for(int dir_z : {-1, 0, 1})
            {
                if (dir_x == 0 && dir_y == 0 && dir_z == 0)
                    continue;

                if (!is_periodic[2] && dir_z) continue;
                if (dir_z < 0 && valid_box.loVect()[2] != domain.loVect()[2]) continue;
                if (dir_z > 0 && valid_box.hiVect()[2] != domain.hiVect()[2]) continue;

                link->add_wrap(diy::Direction{dir_x, dir_y, dir_z});
            }
        }
    }
}

void record_neighbors(int level, int finest_level, const std::vector<int>& gid_offsets, const std::vector<int>& refinements,
        Vector<std::unique_ptr<MultiFab> >& particle_mf, const Box& valid_box,
        const std::array<bool, DIM>& is_periodic, diy::AMRLink* link)
{
    for(int nbr_lev = std::max(0, level - 1); nbr_lev <= std::min(finest_level, level + 1); ++nbr_lev)
    {
        // gotta do this yoga to work around AMReX's static variables
        // TODO: where to get domain? if we have refinements, can
        // calculate it myself from domain of the 0th level in geom_in
        //const Box& nbr_lev_domain; // = particle_mf[level]->boxArray();
        Box nbr_lev_domain; // = particle_mf[level]->boxArray();

        Periodicity periodicity(IntVect(AMREX_D_DECL(nbr_lev_domain.length(0) * is_periodic[0],
                nbr_lev_domain.length(1) * is_periodic[1],
                nbr_lev_domain.length(2) * is_periodic[2])));

        const std::vector<IntVect>& pshifts = periodicity.shiftIntVect();
        // TODO: here we always assume ghosts, get this information somehow
        int ng = 0;
        const BoxArray& ba = particle_mf[nbr_lev]->boxArray();
        // TODO!
//                        int ratio = mesh.RefRatio().at(std::min(lev, nbr_lev));
        int ratio = 2;

        Box gbx = valid_box;
        if (nbr_lev < level)
            gbx.coarsen(ratio);
        else if (nbr_lev > level)
            gbx.refine(ratio);
        gbx.grow(1);

        std::vector<std::pair<int, Box>> isects;

        for(const auto& piv : pshifts)
        {
            ba.intersections(gbx + piv, isects);
            for(const auto& is : isects)
            {
                // is.first is the index of neighbor box
                // ba[is.first] is the neighbor box
                int nbr_gid = gid_offsets.at(nbr_lev) + is.first;
                const Box& nbr_box = ba[is.first];
                Box nbr_ghost_box = grow(nbr_box, ng);
                link->add_neighbor(diy::BlockID{nbr_gid,
                                                -1});        // we don't know the proc, but we'll figure it out later through DynamicAssigner
                link->add_bounds(nbr_lev, refinements.at(nbr_lev), bounds(nbr_box), bounds(nbr_ghost_box));
            }
        }
    } // loop to record neighbors
}


void prepare_master_reader(
        int nblocks,
        diy::Master& master_reader,
        diy::MemoryBuffer& header,
        diy::DiscreteBounds& domain_diy,
        Vector<MultiFab*>& new_state,
        Vector<std::unique_ptr<MultiFab> >& particle_mf,
        int finest_level,
        const amrex::Geometry geom_in)
{
    std::vector<std::string> all_var_names { "particle_mass_density", "density", "xmom", "ymom", "zmom" };

    bool debug = false;
    const int n_levels = finest_level + 1;

    std::array<bool, DIM> is_periodic;
    for(int i = 0; i < DIM; ++i) {
        is_periodic[i] = geom_in.isPeriodic(i);
    }

    const Box& domain = geom_in.Domain();

    for(int i = 0; i < DIM; ++i)
    {
        domain_diy.min[i] = domain.loVect()[i];
        domain_diy.max[i] = domain.hiVect()[i];
    }


    std::vector<int> gid_offsets = {0};
    std::vector<int> refinements = {1};
    nblocks = 0;

    // iterate over all levels to collect refinements and box array sizes (=gid offsets)
    for(int level = 0; level < n_levels; ++level)
    {
        const MultiFab& mf = *new_state[level];

        BoxArray ba = mf.boxArray();
        nblocks += ba.size();
        gid_offsets.push_back(nblocks);

        // TODO: refinements?
        //refinements.push_back(refinements.back() * plotfile.refRatio(level));
    }

    Real* fab_ptr_copy{nullptr};

    std::map<int, Real*> gid_to_fab;
    std::map<int, std::vector<Real*>> gid_to_extra_pointers;

    for(int level = 0; level < n_levels; ++level)
    {
        const MultiFab& mf = *particle_mf[level];
        const BoxArray ba = mf.boxArray();
        BoxArray ba_finer;
        if (level < finest_level)
        {
            // TODO: this will load actual data; we only boxes from finer level
            const MultiFab& mf_finer = *particle_mf[level + 1];
            ba_finer = mf_finer.boxArray();
        }

        // false is for no tiling in MFIter; we want boxes exactly as they are in plotfile
        for(MFIter mfi(mf, false); mfi.isValid(); ++mfi)
        {
            const FArrayBox& dm_fab = mf[mfi];

            const Box& valid_box = mfi.validbox();
            Block::Shape valid_shape;
            for(size_t i = 0; i < DIM; ++i) {
                valid_shape[i] = valid_box.bigEnd()[i] - valid_box.smallEnd()[i] + 1;
            }

            // This is the Box on which the FArrayBox is defined.
            // Note that "abox" includes ghost cells (if there are any),
            // and is thus larger than or equal to "box".
            Box abox = dm_fab.box();

            int gid = gid_offsets[level] + mfi.index();

            std::vector<std::pair<int, Box>> isects;
            diy::AMRLink* link = new diy::AMRLink(3, level, refinements[level], bounds(valid_box), bounds(abox));

            // init fab; dm_fab contains only dark matter density
            Real* fab_ptr = const_cast<Real*>(dm_fab.dataPtr(0));
            long long int fab_size = valid_box.numPts();
            if (valid_box.numPts() != valid_shape[0] * valid_shape[1] * valid_shape[2]) {
                throw std::runtime_error("shape mismatch in valid_box");
            }

            // allocate memory for all fields that we store in FabBlock
            // actual copying for next fields will happen later
            fab_ptr_copy = new Real[fab_size];
            std::vector<Real*> extra_pointers;
            for(int i = 0; i < all_var_names.size(); ++i)
            {
                Real* extra_ptr_copy = new Real[fab_size];
                extra_pointers.push_back(extra_ptr_copy);
            }

            gid_to_fab[gid] = fab_ptr_copy;
            gid_to_extra_pointers[gid] = extra_pointers;

            // copy dark matter density to FabBlock
            memcpy(fab_ptr_copy, fab_ptr, sizeof(Real) * fab_size);
            // copy dark matter density to extra_data
            memcpy(extra_pointers[0], fab_ptr, sizeof(Real) * fab_size);

            master_reader.add(gid, new FabBlockR(fab_ptr_copy, all_var_names, extra_pointers, valid_shape), link);

            record_wrap(domain, valid_box, is_periodic, link);

            record_neighbors(level, finest_level, gid_offsets, refinements, particle_mf, valid_box, is_periodic, link);

            {
                /*
                Real* block_extra_ptr = gid_to_extra_pointers.at(gid).at(var_idx);
                Real* block_fab_ptr = gid_to_fab.at(gid);
                bool add_to_fab = var_idx < n_mt_vars;
                for(int i = 0; i < fab_size; ++i)
                {
                    if (add_to_fab)
                    {
                        block_fab_ptr[i] += fab_ptr[i];
                    }
                    block_extra_ptr[i] = fab_ptr[i];
                }
                */
            }
        } // loop over tiles
    } // loop over levels

    // fill dynamic assigner and fix links
    diy::DynamicAssigner assigner(master_reader.communicator(), master_reader.communicator().size(), nblocks);
    diy::fix_links(master_reader, assigner);

    master_reader.foreach([debug](Block* b, const diy::Master::ProxyWithLink& cp) {
                auto* l = static_cast<diy::AMRLink*>(cp.link());
                auto receivers = link_unique(l, cp.gid());
            }
    );
}

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

void Nyx::runReeberAnalysis(Vector<MultiFab*>& new_state,
        Vector<std::unique_ptr<MultiFab> >& particle_mf,
        const amrex::Geometry geom_in,
        int n_step,
        bool do_analysis,
        std::vector<Halo>& reeber_halos)
{
    amrex::Print() << "ACHTUNG REEBER" << std::endl;

    if(verbose) {
        amrex::Print()<<"Running Reeber anaylsis"<<std::endl;
    }

    if (!do_analysis)
        return;

    int finest_level = parent->finestLevel();

    for (int lev = 0; lev <= parent->finestLevel(); lev++)
    {
        bool do_tiling = false; // TilingIfNotGPU()
        for ( MFIter mfi(*(new_state[lev]), do_tiling); mfi.isValid(); ++mfi )
        {
            const Box& tbx = mfi.tilebox();

            Array4<const Real> state = new_state[lev]->array(mfi);

            Array4<Real> particle = particle_mf[lev]->array(mfi);

            /*
            FArrayBox& fab_state = new_state[lev][mfi];
            FArrayBox& fab_particle = particle_mf[lev][mfi];
            */

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


} // runReeberAnalysis
