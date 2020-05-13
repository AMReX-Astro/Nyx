#include <vector>
#include <string>

#include "Nyx.H"
#include "reeber.H"

#define REEBER_EXTRA_INTEGRAL
#define REEBER_COMPUTE_GAS_VELOCITIES

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
using AmrVertexId = reeber::AmrVertexId;
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

using Grid = Block::Grid;
using GridRef = Block::GridRef;
using Shape = Block::Shape;

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
    for(int i = 0; i < 3; ++i) {
        bounds.min[i] = box.loVect()[i];
        bounds.max[i] = box.hiVect()[i];
    }
    return bounds;
}

void set_wrap(const Box& domain, const Box& valid_box, const std::array<bool, AMREX_SPACEDIM>& is_periodic, diy::AMRLink* link)
{
    for(int dir_x : {-1, 0, 1}) {
        if (!is_periodic[0] && dir_x) continue;
        if (dir_x < 0 && valid_box.loVect()[0] != domain.loVect()[0]) continue;
        if (dir_x > 0 && valid_box.hiVect()[0] != domain.hiVect()[0]) continue;

        for(int dir_y : {-1, 0, 1}) {
            if (!is_periodic[1] && dir_y) continue;
            if (dir_y < 0 && valid_box.loVect()[1] != domain.loVect()[1]) continue;
            if (dir_y > 0 && valid_box.hiVect()[1] != domain.hiVect()[1]) continue;
            for(int dir_z : {-1, 0, 1}) {
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
    for(int nbr_lev = std::max(0, level - 1); nbr_lev <= std::min(finest_level, level + 1); ++nbr_lev) {
        // gotta do this yoga to work around AMReX's static variables

        Box nbr_lev_domain = domain;
        nbr_lev_domain.refine(refinements[level]);

        Periodicity periodicity(IntVect(AMREX_D_DECL(nbr_lev_domain.length(0) * is_periodic[0],
                nbr_lev_domain.length(1) * is_periodic[1],
                nbr_lev_domain.length(2) * is_periodic[2])));

        const std::vector<IntVect>& pshifts = periodicity.shiftIntVect();
        int ng = 0;
        const BoxArray& ba = particle_mf[nbr_lev]->boxArray();
        // TODO: check!
        int ratio;
        if (nbr_lev < level) {
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

        for(const auto& piv : pshifts) {
            ba.intersections(gbx + piv, isects);
            for(const auto& is : isects) {
                // is.first is the index of neighbor box
                // ba[is.first] is the neighbor box
                int nbr_gid = gid_offsets.at(nbr_lev) + is.first;
                const Box& nbr_box = ba[is.first];
                Box nbr_ghost_box = grow(nbr_box, ng);
                link->add_neighbor(diy::BlockID{nbr_gid,
                                                -1});        // we don't know the proc, but we'll figure it out later through DynamicAssigner
                // we ignore ghost data, hence nbr_box in last two parameter
                link->add_bounds(nbr_lev, refinements.at(nbr_lev), bounds(nbr_box), bounds(nbr_box));
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
        const amrex::Geometry geom_in,
        std::vector<std::unique_ptr<Real[]>>& pointers_to_copied_data)
{
    std::vector<std::string> new_state_vars { "density", "xmom", "ymom", "zmom" };
    const size_t density_var_idx = 0;

    std::vector<std::string> all_vars { "particle_mass_density", "density", "xmom", "ymom", "zmom" };

    bool debug = false;
    const int n_levels = finest_level + 1;

    std::array<bool, AMREX_SPACEDIM> is_periodic;
    for(int i = 0; i < AMREX_SPACEDIM; ++i) {
        is_periodic[i] = geom_in.isPeriodic(i);
    }

    const Box& domain = geom_in.Domain();

    for(int i = 0; i < AMREX_SPACEDIM; ++i) {
        diy_domain.min[i] = domain.loVect()[i];
        diy_domain.max[i] = domain.hiVect()[i];
    }

    std::vector<int> gid_offsets = {0};
    std::vector<int> refinements = {1};
    int nblocks = 0;

    // iterate over all levels to collect refinements and box array sizes (=gid offsets)
    for(int level = 0; level < n_levels; ++level) {
        const MultiFab& mf = *new_state[level];

        BoxArray ba = mf.boxArray();
        nblocks += ba.size();
        gid_offsets.push_back(nblocks);

        if (level > 0) {
            // we accumulate in refinements ratio to the coarsest level
            // assuming uniform refinement ratio in all dimensions
            // level_refinements contains fine_ratio
            refinements.push_back(refinements.back() * level_refinements[level][0]);
        }
    }

    //fmt::print("REFINEMENTS = {}\n", container_to_string(refinements));

    std::map<int, Real*> gid_to_fab;
    std::map<int, long long int> gid_to_fab_size;
    std::map<int, std::vector<Real*>> gid_to_extra_pointers;

    for(int level = 0; level < n_levels; ++level) {
        const MultiFab& dm_mf = *particle_mf[level];
        const BoxArray ba = dm_mf.boxArray();

        // false is for no tiling in MFIter; we want boxes exactly as they are in plotfile
        for(MFIter mfi(dm_mf, false); mfi.isValid(); ++mfi) {
            const FArrayBox& dm_fab = dm_mf[mfi];
            int gid = gid_offsets[level] + mfi.index();


            // core stuff
            const Box& core_box = mfi.validbox();
            Shape core_shape;
            for(size_t i = 0; i < AMREX_SPACEDIM; ++i) {
                core_shape[i] = core_box.bigEnd()[i] - core_box.smallEnd()[i] + 1;
            }
            long long int core_fab_size = core_box.numPts();
            gid_to_fab_size[gid] = core_fab_size;

            pointers_to_copied_data.emplace_back(new Real[core_fab_size]);
            Real* core_fab_ptr = pointers_to_copied_data.back().get();

            GridRef core_grid_ref(core_fab_ptr, core_shape, false);

            // ghost stuff
            Box ghost_box = dm_fab.box();
            Shape ghost_shape;
            for(size_t i = 0; i < AMREX_SPACEDIM; ++i) {
                ghost_shape[i] = ghost_box.bigEnd()[i] - ghost_box.smallEnd()[i] + 1;
            }
            // this grid ref points directly to dm_fab data
            // index 0, since dm_fab contains only dark matter density
            GridRef ghost_grid_ref(const_cast<Real*>(dm_fab.dataPtr(0)), ghost_shape, /* c_order = */ false);

            Shape core_ghost_adjustment;
            for(size_t i = 0; i < AMREX_SPACEDIM; ++i) {
                core_ghost_adjustment[i] = core_box.smallEnd()[i] - ghost_box.smallEnd()[i];
                if (core_ghost_adjustment[i] < 0) { throw std::runtime_error("ghost box must be larger than core box"); }
            }

            // copy particle data in the core to core_grid_ref
            // reason: dm_fab has 2 ghosts, new_state has 1, and we don't need
            // ghosts
            diy::for_each(core_shape, [&core_grid_ref, &ghost_grid_ref, core_ghost_adjustment](const Vertex& v) {
                    core_grid_ref(v) = ghost_grid_ref(v + core_ghost_adjustment);
                    });

            std::vector<std::pair<int, Box>> isects;
            diy::AMRLink* link = new diy::AMRLink(3, level, refinements[level], bounds(core_box), bounds(core_box));


            //if (debug) fmt::print(std::cerr, "Added box, particle_mf ghost_box = {}, core_box = {}, core_fab_size = {}\n", ghost_box, core_box, core_fab_size);

            // allocate memory for all fields that we store in FabBlock
            // actual copying for next fields will happen later
            std::vector<Real*> extra_pointers;

            // reserve memory for dm_density
            pointers_to_copied_data.emplace_back(new Real[core_fab_size]);
            Real* extra_ptr_copy = pointers_to_copied_data.back().get();
            extra_pointers.push_back(extra_ptr_copy);

            // reserve memory for variables in new_state
            for(int i = 0; i < new_state_vars.size(); ++i) {
                Real* extra_ptr_copy = new Real[core_fab_size];
                extra_pointers.push_back(extra_ptr_copy);
            }

            gid_to_fab[gid] = core_fab_ptr;
            gid_to_extra_pointers[gid] = extra_pointers;

            // copy dark matter density to extra_data
            memcpy(extra_pointers[0], core_fab_ptr, sizeof(Real) * core_fab_size);

            if (all_vars.size() != extra_pointers.size())
                throw std::runtime_error("all_vars.size() != extra_pointers.size()");

            master_reader.add(gid, new FabBlockR(core_fab_ptr, all_vars, extra_pointers, core_shape), link);

            set_wrap(domain, core_box, is_periodic, link);

            set_neighbors(level, finest_level, gid_offsets, refinements, domain, particle_mf, core_box, is_periodic, link);
        } // loop over tiles
    } // loop over levels for dark matter

    // now FabBlocks are created, it remains to add gas density
    // and copy velocities to FabBlocks
    for(int level = 0; level < n_levels; ++level) {
        const MultiFab& state_mf = *new_state[level];
        // false is for no tiling in MFIter; we want boxes exactly as they are
        for(MFIter mfi(state_mf, false); mfi.isValid(); ++mfi) {
            const FArrayBox& state_fab = state_mf[mfi];

            int gid = gid_offsets[level] + mfi.index();
            long long int fab_size = gid_to_fab_size.at(gid);

            //fmt::print(std::cerr, "new_state_vars, level = {}, gid = {}, fab_size = {}, mfi.index = {}, gid[offsets[level] = {}\n",
                    //level, gid, fab_size, mfi.index(), gid_offsets[level]);

            // core box and shape
            const Box& core_box = mfi.validbox();
            if (core_box.numPts() != fab_size) { throw std::runtime_error("shape mismatch state_mf != particle_mf"); }
            Shape core_shape;
            for(size_t i = 0; i < AMREX_SPACEDIM; ++i) { core_shape[i] = core_box.length(i); }

            // ghost box and shape
            Box ghost_box = state_fab.box();
            Shape ghost_shape;
            for(size_t i = 0; i < AMREX_SPACEDIM; ++i) { ghost_shape[i] = ghost_box.length(i); }

            Shape core_ghost_adjustment;
            for(size_t i = 0; i < AMREX_SPACEDIM; ++i) { core_ghost_adjustment[i] = core_box.smallEnd()[i] - ghost_box.smallEnd()[i]; }

            //if (debug) fmt::print(std::cerr, "state_fab  ghost_box = {}, core_box = {}, state_mf.contains_nan = {} (grow = 0)\n",
            //        ghost_box, core_box, state_mf.contains_nan(0, new_state_vars.size(), 0, false));

            size_t gas_density_idx = 0;

            for(int var_idx = 0; var_idx < new_state_vars.size(); ++var_idx) {

                bool add_to_fab = (var_idx == 0);

                // this grid ref points directly to dm_fab data
                GridRef ghost_grid_ref(const_cast<Real*>(state_fab.dataPtr(var_idx)), ghost_shape, /* c_order = */ false);

                Real* block_extra_ptr = gid_to_extra_pointers.at(gid).at(var_idx);
                GridRef core_grid_ref(block_extra_ptr, core_shape, false);

                if (add_to_fab) {
                    // density - copy to extra_data as is, add to main grid

                   // copy core data to block_extra_ptr
                    diy::for_each(core_shape, [&core_grid_ref, &ghost_grid_ref, core_ghost_adjustment](const Vertex& v) {
                            core_grid_ref(v) = ghost_grid_ref(v + core_ghost_adjustment);
                            });

                    //if (debug) fmt::print(std::cerr, "copied variable {} to extra_pointers, add_to_fab = {}\n",
                    //        new_state_vars[var_idx], add_to_fab);

                    Real* block_fab_ptr = gid_to_fab.at(gid);
                    GridRef core_grid_ref(block_extra_ptr, core_shape, false);
                    // add core density data to block_fab_ptr
                    diy::for_each(core_shape, [&core_grid_ref, &ghost_grid_ref, core_ghost_adjustment](const Vertex& v) {
                            core_grid_ref(v)  += ghost_grid_ref(v + core_ghost_adjustment);
                            });

                } else {
                    // momenta - divide by density

                    GridRef density_grid_ref(const_cast<Real*>(state_fab.dataPtr(density_var_idx)), ghost_shape, /* c_order = */ false);

                    diy::for_each(core_shape, [&core_grid_ref, &ghost_grid_ref, &density_grid_ref, core_ghost_adjustment](const Vertex& v) {
                            auto ghost_v = v + core_ghost_adjustment;
                            core_grid_ref(v) = ghost_grid_ref(ghost_v) / density_grid_ref(ghost_v);
                            });

                    //if (debug) fmt::print(std::cerr, "copied {} to extra_pointers, divided by density\n",
                    //        new_state_vars[var_idx], add_to_fab);
                } // if - density or momenta
            } // loop over vars
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
                                diy::DiscreteBounds diy_domain,
                                Real absolute_rho,
                                bool negate,
                                Real min_halo_n_cells)
{
    bool debug = false; //world.rank() == 0;
    std::string prefix = "./DIY.XXXXXX";
    int in_memory = -1;
    std::string log_level = "info";
    diy::FileStorage storage(prefix);

    std::vector<Halo> result;

    //if (debug) fmt::print(std::cerr, "Master reader - started\n");

    // copy FabBlocks to FabComponentBlocks
    // in FabTmtConstructor mask will be set and local trees will be computed
    // FabBlock can be safely discarded afterwards
    diy::Master master(world, threads, in_memory, &Block::create, &Block::destroy, &storage, &Block::save,
            &Block::load);

    //if (debug) fmt::print(std::cerr, "Master reader created\n");

    Real cell_volume = geom_in.ProbDomain().volume() / geom_in.Domain().numPts();


    //if (debug) fmt::print(std::cerr, "Cell volume = {}\n", cell_volume);

    master_reader.foreach(
            [&master, debug, diy_domain, absolute_rho, negate, cell_volume](FabBlockR* b, const diy::Master::ProxyWithLink& cp) {
                auto* l = static_cast<AMRLink*>(cp.link());
                AMRLink* new_link = new AMRLink(*l);

                // prepare neighbor box info to save in MaskedBox
                // TODO: refinement vector
                int local_ref = l->refinement()[0];
                int local_lev = l->level();

                //if (debug) fmt::print(std::cerr, "adding block, bounds = {}, core = {}, gid = {}\n", l->bounds(), l->core(), cp.gid());

                master.add(cp.gid(),
                        new Block(b->fab, b->extra_names_, b->extra_fabs_, local_ref, local_lev, diy_domain,
                                l->bounds(),
                                l->core(), cp.gid(),
                                new_link, absolute_rho, negate, /*absolute = */ true, cell_volume),
                        new_link);

                //if (debug) fmt::print(std::cerr, "Added block gid = {}\n", cp.gid());
            });

    //if (debug) fmt::print(std::cerr, "Master populated\n");

    int global_n_undone = 1;

    master.foreach(&send_edges_to_neighbors_cc<Real, AMREX_SPACEDIM>);
    master.exchange();
    master.foreach(&delete_low_edges_cc<Real, AMREX_SPACEDIM>);

    //if (debug) fmt::print(std::cerr, "Low edges deleted\n");

    int rounds = 0;
    while(global_n_undone) {
        rounds++;

        master.foreach(&amr_cc_send<Real, AMREX_SPACEDIM>);
        master.exchange();
        master.foreach(&amr_cc_receive<Real, AMREX_SPACEDIM>);
        master.exchange();

        global_n_undone = master.proxy(master.loaded_block()).read<int>();

        //if (debug) fmt::print(std::cerr, "global_n_undone = {}\n", global_n_undone);

    }

    master.foreach([](Block* b, const diy::Master::ProxyWithLink& cp) {
        b->compute_final_connected_components();
        b->compute_local_integral();
    });

    //if (debug) fmt::print(std::cerr, "Local integrals computed");

    amrex::Box domain = geom_in.Domain();
    RealBox prob_domain = geom_in.ProbDomain();

    master.foreach(
            [&result, domain, prob_domain, min_halo_n_cells](
                    Block* b,
                    const diy::Master::ProxyWithLink& cp) {

                diy::Point<int, AMREX_SPACEDIM> domain_shape;
                Real domain_volume = domain.numPts();
                for(int i = 0; i < AMREX_SPACEDIM; ++i) {
                    domain_shape[i] = domain.length(i);
                }

                diy::GridRef<void*, 3> domain_box(nullptr, domain_shape, /* c_order = */ false);

                for(auto& root_values_pair : b->local_integral_) {
                    AmrVertexId root = root_values_pair.first;
                    if (root.gid != b->gid)
                        continue;

                    auto& values = root_values_pair.second;

                    Real halo_n_cells = values.at("n_cells");

                    if (halo_n_cells < min_halo_n_cells)
                        continue;

                    Halo new_halo;

                    const Real cell_volume = b->cell_volume_;

                    new_halo.n_cells = halo_n_cells;
                    new_halo.volume = (Real)halo_n_cells * cell_volume;
                    new_halo.gas_mass = cell_volume * values["density"];
                    new_halo.dm_mass = cell_volume * values["particle_mass_density"];
                    new_halo.total_mass = new_halo.gas_mass + new_halo.dm_mass;

                    // division by density was done in prepare_master_reader,
                    // it remains to multiply with volume to get velocities
                    new_halo.gas_vel_x = cell_volume * values["xmom"];
                    new_halo.gas_vel_y = cell_volume * values["ymom"];
                    new_halo.gas_vel_z = cell_volume * values["zmom"];

                    auto halo_root_position = b->local_.global_position(root);

                    int ref = b->refinement();
                    amrex::IntVect ref_ratio (ref);
                    amrex::Box domain_at_level = domain;
                    domain_at_level.refine(ref_ratio);

                    new_halo.id = domain_box.index(halo_root_position);

                    for(int i = 0; i < AMREX_SPACEDIM; ++i) {
                        new_halo.position[i] = halo_root_position[i];

                        new_halo.real_position[i] = prob_domain.lo(i) +
                            prob_domain.length(i) * (Real)(halo_root_position[i]) / domain_at_level.length(i);
                    }
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
    if (!do_analysis)
        return;

    verbose = true;

    BL_PROFILE("Nyx::runReeberAnalysis()");

    if(verbose) {
        amrex::Print() <<   "Running Reeber anaylsis" << std::endl;
    }

    // store pointers to all dynamically allocated arrays, so that
    // data will be freed automatically after exiting compute_halos
    std::vector<std::unique_ptr<Real[]>> pointers_to_copied_data;

    diy::mpi::communicator world = ParallelDescriptor::Communicator();

    int threads = 1;
    std::string prefix = "./DIY.XXXXXX";
    int in_memory = -1;
    diy::FileStorage storage(prefix);

    diy::MemoryBuffer header;

    BL_PROFILE_VAR("Nyx::runReeberAnalysis()::master_reader",master_reader_var);

    diy::Master master_reader(world, threads, in_memory, &FabBlockR::create, &FabBlockR::destroy,
            &storage, &FabBlockR::save, &FabBlockR::load);

    BL_PROFILE_VAR_STOP(master_reader_var);

    diy::DiscreteBounds diy_domain(3);

    // TODO: take rho, min_halo_n_cells as parameters
    Real min_halo_n_cells = 10;
    Real rho = 81.66;

    Real absolute_rho = (Nyx::average_dm_density + Nyx::average_gas_density) * rho;
    bool negate = true;  // sweep superlevel sets, highest density = root
    int finest_level = parent->finestLevel();

    //if (verbose and world.rank() == 0) fmt::print(std::cerr, "fmt prepare_master_reader called, finest_level = {}, rho = {}\n", finest_level, absolute_rho);

    BL_PROFILE_VAR("Nyx::runReeberAnalysis()::prepare_master_reader",prepare_master_reader_var);

    prepare_master_reader(master_reader, header, diy_domain, new_state, particle_mf,
                          finest_level, level_refinements, geom_in, pointers_to_copied_data);

    BL_PROFILE_VAR_STOP(prepare_master_reader_var);

    //if (verbose and world.rank() == 0) fmt::print(std::cerr, "prepare_master_reader finished\n");

    BL_PROFILE_VAR("Nyx::runReeberAnalysis()::compute_halos",compute_halos_var);

    reeber_halos = compute_halos(world, master_reader, geom_in, threads, diy_domain, absolute_rho, negate, min_halo_n_cells);
    if (verbose and world.rank() == 0) fmt::print(std::cerr, "compute_halos finished, result.size = {}\n", reeber_halos.size());

    BL_PROFILE_VAR_STOP(compute_halos_var);

} // runReeberAnalysis
