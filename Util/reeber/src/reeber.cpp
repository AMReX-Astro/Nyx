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

#include <amr-connected-components-complex.h>


using namespace amrex;

// block-independent types
using AMRLink = diy::AMRLink;

using Bounds = diy::DiscreteBounds;
using AmrVertexId = r::AmrVertexId;
using AmrEdge = reeber::AmrEdge;

using FabBlockR = FabBlock<Real, AMREX_SPACEDIM>;

using Block = FabComponentBlock<Real, AMREX_SPACEDIM>;
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

void set_wrap(const Box& domain, const Box& valid_box, const std::array<bool, AMREX_SPACEDIM>& is_periodic, diy::AMRLink* link)
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

void set_neighbors(int level, int finest_level, const std::vector<int>& gid_offsets, const std::vector<int>& refinements,
        const Box& domain, Vector<std::unique_ptr<MultiFab> >& particle_mf, const Box& valid_box,
        const std::array<bool, AMREX_SPACEDIM>& is_periodic, diy::AMRLink* link)
{
    for(int nbr_lev = std::max(0, level - 1); nbr_lev <= std::min(finest_level, level + 1); ++nbr_lev)
    {
        // gotta do this yoga to work around AMReX's static variables

        Box nbr_lev_domain = domain;
        nbr_lev_domain.refine(refinements[level]);

        Periodicity periodicity(IntVect(AMREX_D_DECL(nbr_lev_domain.length(0) * is_periodic[0],
                nbr_lev_domain.length(1) * is_periodic[1],
                nbr_lev_domain.length(2) * is_periodic[2])));

        const std::vector<IntVect>& pshifts = periodicity.shiftIntVect();
        // TODO: here we always assume no ghosts, get this information somehow
        int ng = 0;
        const BoxArray& ba = particle_mf[nbr_lev]->boxArray();
        // TODO: check!
        int ratio;
        if (nbr_lev < level)
        {
            ratio = refinements[level] / refinements[nbr_lev];
        } else {
            ratio = refinements[nbr_lev] / refinements[level];
        }

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
    } // loop to set neighbors
}


void prepare_master_reader(
        diy::Master& master_reader,
        diy::MemoryBuffer& header,
        diy::DiscreteBounds& diy_domain,
        Vector<MultiFab*>& new_state,
        Vector<std::unique_ptr<MultiFab> >& particle_mf,
        int finest_level,
        const Vector<IntVect>& level_refinements,
        const amrex::Geometry geom_in)
{
    std::vector<std::string> new_state_vars { "density", "xmom", "ymom", "zmom" };
    std::vector<std::string> all_vars { "particle_mass_density", "density", "xmom", "ymom", "zmom" };

    bool debug = false;
    const int n_levels = finest_level + 1;

    std::array<bool, AMREX_SPACEDIM> is_periodic;
    for(int i = 0; i < AMREX_SPACEDIM; ++i) {
        is_periodic[i] = geom_in.isPeriodic(i);
    }

    const Box& domain = geom_in.Domain();

    for(int i = 0; i < AMREX_SPACEDIM; ++i)
    {
        diy_domain.min[i] = domain.loVect()[i];
        diy_domain.max[i] = domain.hiVect()[i];
    }

    std::vector<int> gid_offsets = {0};
    std::vector<int> refinements = {1};
    int nblocks = 0;

    // iterate over all levels to collect refinements and box array sizes (=gid offsets)
    for(int level = 0; level < n_levels; ++level)
    {
        const MultiFab& mf = *new_state[level];

        BoxArray ba = mf.boxArray();
        nblocks += ba.size();
        gid_offsets.push_back(nblocks);

        // we accumulate in refinements ratio to the coarsest level
        // assuming uniform refinement ratio in all dimensions
        // level_refinements contains fine_ratio
        refinements.push_back(refinements.back() * level_refinements[level][0]);
    }

    Real* fab_ptr_copy{nullptr};

    std::map<int, Real*> gid_to_fab;
    std::map<int, long long int> gid_to_fab_size;
    std::map<int, std::vector<Real*>> gid_to_extra_pointers;

    for(int level = 0; level < n_levels; ++level)
    {
        const MultiFab& dm_mf = *particle_mf[level];
        const BoxArray ba = dm_mf.boxArray();

        // false is for no tiling in MFIter; we want boxes exactly as they are in plotfile
        for(MFIter mfi(dm_mf, false); mfi.isValid(); ++mfi)
        {
            const FArrayBox& dm_fab = dm_mf[mfi];

            const Box& valid_box = mfi.validbox();
            Block::Shape valid_shape;
            for(size_t i = 0; i < AMREX_SPACEDIM; ++i) {
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
            gid_to_fab_size[gid] = fab_size;
            if (valid_box.numPts() != valid_shape[0] * valid_shape[1] * valid_shape[2]) {
                throw std::runtime_error("shape mismatch in valid_box");
            }

            // allocate memory for all fields that we store in FabBlock
            // actual copying for next fields will happen later
            fab_ptr_copy = new Real[fab_size];
            std::vector<Real*> extra_pointers;

            // reserve memory for dm_density
            Real* extra_ptr_copy = new Real[fab_size];
            extra_pointers.push_back(extra_ptr_copy);

            // reserve memory for variables in new_state
            for(int i = 0; i < new_state_vars.size(); ++i)
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

            if (all_vars.size() != extra_pointers.size())
                throw std::runtime_error("all_vars.size() != extra_pointers.size()");

            master_reader.add(gid, new FabBlockR(fab_ptr_copy, all_vars, extra_pointers, valid_shape), link);

            set_wrap(domain, valid_box, is_periodic, link);

            set_neighbors(level, finest_level, gid_offsets, refinements, domain, particle_mf, valid_box, is_periodic, link);
        } // loop over tiles
    } // loop over levels for dark matter

    for(int level = 0; level < n_levels; ++level)
    {
        const MultiFab& state_mf = *new_state[level];
        // false is for no tiling in MFIter; we want boxes exactly as they are
        for(MFIter mfi(state_mf, false); mfi.isValid(); ++mfi)
        {
            const FArrayBox& state_fab = state_mf[mfi];

            int gid = gid_offsets[level] + mfi.index();
            long long int fab_size = gid_to_fab_size[gid];
            for(int var_idx = 0; var_idx < new_state_vars.size(); ++var_idx)
            {
                Real* fab_ptr = const_cast<Real*>(state_fab.dataPtr(var_idx));
                Real* block_extra_ptr = gid_to_extra_pointers.at(gid).at(var_idx);
                Real* block_fab_ptr = gid_to_fab.at(gid);
                bool add_to_fab = new_state_vars[var_idx] == "density";
                // momentum will be divied by density in FabComponentBlock ctor
                // later, here we just copy it
                for(int i = 0; i < fab_size; ++i)
                {
                    if (add_to_fab)
                    {
                        block_fab_ptr[i] += fab_ptr[i];
                    }
                    block_extra_ptr[i] = fab_ptr[i];
                }
            }
        } // loop over tiles
    } // loop over levels for new state


    // fill dynamic assigner and fix links
    diy::DynamicAssigner assigner(master_reader.communicator(), master_reader.communicator().size(), nblocks);
    diy::fix_links(master_reader, assigner);

    master_reader.foreach([](FabBlockR* b, const diy::Master::ProxyWithLink& cp) {
                auto* l = static_cast<diy::AMRLink*>(cp.link());
                auto receivers = link_unique(l, cp.gid());
            }
    );
}


std::vector<Halo> compute_halos(diy::mpi::communicator& world,
                                diy::Master& master_reader,
                                const amrex::Geometry geom_in,
                                int threads,
                                diy::DiscreteBounds domain,
                                Real absolute_rho,
                                bool negate,
                                Real min_halo_volume)
{
    std::string prefix = "./DIY.XXXXXX";
    int in_memory = -1;
    std::string log_level = "info";
    diy::FileStorage storage(prefix);

    std::vector<Halo> result;

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

                // after master initialized its blocks,
                // data in b->fab is no longer needed;
                // extra_fabs will be released in FabComponentBlock::init
                delete [] b->fab.data();

            });

    int global_n_undone = 1;

    master.foreach(&send_edges_to_neighbors_cc<Real, AMREX_SPACEDIM>);
    master.exchange();
    master.foreach(&delete_low_edges_cc<Real, AMREX_SPACEDIM>);

    int rounds = 0;
    while(global_n_undone)
    {
        rounds++;

        master.foreach(&amr_cc_send<Real, AMREX_SPACEDIM>);
        master.exchange();
        master.foreach(&amr_cc_receive<Real, AMREX_SPACEDIM>);
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
    bool has_momentum = true;

    RealBox prob_domain = geom_in.ProbDomain();

    master.foreach(
            [&result, domain, prob_domain, min_halo_volume,
                    has_density, has_particle_mass_density,
                    has_momentum](
                    Block* b,
                    const diy::Master::ProxyWithLink& cp) {

                diy::Point<int, AMREX_SPACEDIM> domain_shape;
                Real domain_volume = 1;
                for(int i = 0; i < AMREX_SPACEDIM; ++i)
                {
                    domain_shape[i] = domain.max[i] - domain.min[i] + 1;
                    domain_volume *= domain_shape[i];
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

                    Halo new_halo;

                    // TODO: position in multi-level case
                    new_halo.int_volume = halo_volume;
                    new_halo.volume = (Real)halo_volume * prob_domain.volume() / domain_volume;

                    auto halo_root_position = b->local_.global_position(root);

                    for(int i = 0; i < AMREX_SPACEDIM; ++i) {
                        new_halo.position[i] = halo_root_position[i];

                        new_halo.real_position[i] = prob_domain.lo(i) +
                            (prob_domain.hi(i) - prob_domain.lo(i)) * (Real)(halo_root_position[i]) / (Real)(domain_shape[i]);
                    }

                    new_halo.id = domain_box.index(halo_root_position);
                    new_halo.gas_mass = m_gas;
                    new_halo.dm_mass = m_particles;
                    new_halo.total_mass = m_total;

                    result.push_back(new_halo);
                }
            });

    return result;
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
        const Vector<IntVect>& level_refinements,
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


    diy::mpi::communicator world = ParallelDescriptor::Communicator();

    int threads = 1;
    std::string prefix = "./DIY.XXXXXX";
    int in_memory = -1;
    diy::FileStorage storage(prefix);

    diy::MemoryBuffer header;

    diy::Master master_reader(world, threads, in_memory, &FabBlockR::create, &FabBlockR::destroy,
            &storage, &FabBlockR::save, &FabBlockR::load);

    diy::DiscreteBounds diy_domain;

    // TODO: get threshold
    Real absolute_rho = 81.66;
    bool negate = true;  // sweep superlevel sets, highest density = root
    // TODO: take as parameter
    Real min_halo_volume = 10;

    int finest_level = parent->finestLevel();

    prepare_master_reader(master_reader, header, diy_domain, new_state, particle_mf,
                          finest_level, level_refinements, geom_in);

    reeber_halos = compute_halos(world, master_reader, geom_in, threads, diy_domain, absolute_rho, negate, min_halo_volume);
} // runReeberAnalysis
