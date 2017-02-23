
#ifndef _Nyx_H_
#define _Nyx_H_

#include <BC_TYPES.H>
#include <AmrLevel.H>
#include <ErrorList.H>
#include <FluxRegister.H>

#include "NyxParticleContainer.H"

#include <iostream>

using std::istream;
using std::ostream;

typedef NyxParticleContainer<1+BL_SPACEDIM+3> StellarParticleContainer;
typedef NyxParticleContainer<1+BL_SPACEDIM+3> AGNParticleContainer;

#ifdef NO_HYDRO
enum StateType {
    PhiGrav_Type = 0,
    Gravity_Type,
    NUM_STATE_TYPE
};
#endif

#ifndef NO_HYDRO
enum StateType {
    State_Type = 0,
    DiagEOS_Type,
#ifdef GRAVITY
    PhiGrav_Type,
    Gravity_Type,
#endif
    NUM_STATE_TYPE
};
#endif

//
// AmrLevel-derived class for hyperbolic conservation equations for stellar
// media
//

class Nyx
    :
    public AmrLevel
{
public:
    //
    //Default constructor.  Builds invalid object.
    //
    Nyx();

    //
    //The basic constructor.
    //
    Nyx(Amr& papa, int lev, const Geometry& level_geom, const BoxArray& bl,
        Real time);

    //
    //The destructor.
    //
    virtual ~Nyx();

    //
    //Restart from a checkpoint file.
    //
    virtual void restart(Amr& papa, istream& is, bool b_read_special=false);

    //
    //Call AmrLevel::checkPoint and then add radiation info
    //
    virtual void checkPoint(const std::string& dir, std::ostream& os,
                            VisMF::How how, bool dump_old);

    // A string written as the first item in `write_plot_file()` at level zero.
    // It is so we can distinguish between different types of plot files. For
    // Nyx it has the form: Nyx-Vnnn.
    virtual std::string thePlotFileType() const;

    virtual void setPlotVariables();

    //
    //Write a plotfile to specified directory.
    //
    virtual void writePlotFile(const std::string& dir, ostream& os, VisMF::How how);

    //
    //Write MultiFab as plot file
    //
    void writeMultiFabAsPlotFile(const std::string& pltfile, const MultiFab& mf, std::string componentName);
    //
    //Write all parameters into specified directory.
    //
    static int write_parameters_in_plotfile;
    virtual void write_parameter_file(const std::string& dir);

    //
    // If true then print the warnings from the Fortran routines
    //
    static int print_fortran_warnings;

    //
    //Define data descriptors.
    //
    static void variable_setup();
    static void variable_setup_for_new_comp_procs();
    static void hydro_setup();
    static void no_hydro_setup();

    //
    //Define tagging functions.
    //
    static void error_setup();
    //
    //Cleanup data descriptors at end of run.
    //
    static void variable_cleanup();

    //
    //Initialize grid data at problem start-up.
    //
    virtual void initData();

    //
    //Initialize grid data from a plotfile at problem start-up.
    //
    void init_from_plotfile();

    //
    //Initialize grid data from a plotfile at problem start-up.
    //
    void ReadPlotFile(bool first, const std::string& plot_file_name, bool& rhoe_infile);

    //
    // Read inputs related to comoving coordinates
    //
    static void read_comoving_params();

    //
    // Initial redshift (default = 0)
    //
    static Real initial_z;

    //
    // Final a -- used as stopping criterion if positive
    //
    static Real final_a;

    //
    // Final z -- used as stopping criterion if positive
    //
    static Real final_z;

    //
    // Relative change in a allowed in one timestep
    //
    static Real relative_max_change_a;

    //
    // Old and new times at which "old_a" and "new_a" are defined.
    //
    static Real old_a_time;
    static Real new_a_time;

    //
    // "a" at old_a_time and new_a_time
    //
    static Real old_a;
    static Real new_a;

    //
    // Get the comoving coordinate "a"
    //
    Real get_comoving_a(Real time);
    void integrate_comoving_a(Real time, Real dt);

    //
    // Time step control based on "a" not growing too fast
    //
    void comoving_est_time_step(Real& cur_time, Real& est_dt);

    //
    // Time step control based on "z" not passing one of the specified plot_z_values
    //
    void plot_z_est_time_step(Real& est_dt, bool& dt_changed);

    //
    // Time step control based on "z" not passing one of the specified analysis_z_values
    //
    void analysis_z_est_time_step(Real& est_dt, bool& dt_changed);

    //
    // How to initialize "a" at restart (from checkpoint or plotfile)
    //
    void comoving_a_post_restart(const std::string& restart_file);

#ifdef GRAVITY
    //
    // Initialize the Santa Barbara problem
    //
    void init_santa_barbara(int init_sb_vels);
#endif

    //
    // Initialize from MUSIC
    //
    void initcosmo();
#ifdef GRAVITY
    void icReadAndPrepareFab(std::string mfDirName, int nghost, MultiFab &mf);
#endif

    //
    // Read particle-related inputs
    //
    static void read_particle_params();

    //
    // Read initialization-related inputs
    //
    static void read_init_params();

    //
    // Write particles in checkpoint directories
    //
    void particle_check_point(const std::string& dir);

    //
    // Write particles in plotfile directories
    //
    void particle_plot_file(const std::string& dir);

    //
    // How to initialize at restart
    //
    void particle_post_restart(const std::string& restart_file, bool is_checkpoint = true);

    //
    // Redistribute
    //
    void particle_redistribute(int lbase = 0, bool init = false);

    //
    // Move randomly
    //
    void particle_move_random();

    //
    // Initialize particle locations and velocities (and strengths if relevant)
    //
    virtual void init_particles();

    //
    // Setup virtual particles if necessary
    //
    void setup_virtual_particles();

    //
    // Remove virtual particles if necessary
    //
    void remove_virtual_particles();
    //
    // Setup ghost particles (for finer levels) if necessary
    //
    void setup_ghost_particles();

    //
    // Remove ghost particles (for this level) if necessary
    //
    void remove_ghost_particles();

    //
    // Time step control based on particles
    //
    void particle_est_time_step(Real& est_dt);

    //
    // Derived quantities associated with particles
    //
    MultiFab* particle_derive (const std::string& name, Real time, int ngrow);

    //
    // Default verbosity of Particle class
    //
    static int particle_verbose;

    //
    // Default cfl of particles in Particle class
    //
    static Real particle_cfl;
#ifdef NEUTRINO_PARTICLES
    static Real neutrino_cfl;
#endif

    //
    // Should we write particles into plotfiles?
    //
    static int write_particles_in_plotfile;

    //
    // Shall we write the initial single-level particle density into a multifab
    //   called "ParticleDensity"?
    //
    static int write_particle_density_at_init;

    //
    // Shall we write an ascii file with only the particles in every other cell --
    //   this is a cheap way of creating a "coarser" particle file
    //
    static int write_coarsened_particles;

    //
    //Set time levels of state data.
    //
    virtual void setTimeLevel (Real time, Real dt_old, Real dt_new);

    //
    //Initialize data on this level from another Nyx (during regrid).
    //
    virtual void init(AmrLevel& old);

    //
    // Initialize data on this level after regridding if old level did not
    // previously exist
    //
    virtual void init();

    //
    // Proceed with next timestep?
    //
    virtual int okToContinue();

    //
    // Tell Amr to write a plotfile now
    //
    bool writePlotNow();

    //
    // Tell Nyx to do analysis now
    //
    bool doAnalysisNow();

    //
    // Advance grids at this level in time.
    //
    virtual Real advance(Real time, Real dt, int iteration, int ncycle);

    Real advance_hydro(Real time, Real dt, int iteration, int ncycle);
    Real advance_no_hydro(Real time, Real dt, int iteration, int ncycle);
    Real advance_hydro_plus_particles(Real time, Real dt, int iteration,
                                      int ncycle);
    void just_the_hydro(Real time, Real dt, Real a_old, Real a_new);
    void strang_first_step  (Real time, Real dt,  MultiFab& state, MultiFab&  dstate);
    void strang_second_step (Real time, Real dt,  MultiFab& state, MultiFab&  dstate);

    Real advance_particles_only (Real time, Real dt, int iteration, int ncycle);

    void moveKickDriftExact(Real dt);
    void moveKickExact(Real dt);
    void time_center_source_terms(MultiFab& S_new, MultiFab& ext_src_old,
                                  MultiFab& ext_src_new, Real dt);
    //
    // Estimate time step.
    //
    Real est_time_step(Real dt_old);

    //
    // Compute initial time step.
    //
    Real initial_time_step();

    //
    // Compute initial `dt'.
    //
    virtual void computeInitialDt(int finest_level, int sub_cycle,
                                  Array<int>& n_cycle,
                                  const Array<IntVect>& ref_ratio,
                                  Array<Real>& dt_level, Real stop_time);
    //
    // Compute new `dt'.
    //
    virtual void computeNewDt(int finest_level, int sub_cycle,
                              Array<int>& n_cycle,
                              const Array<IntVect>& ref_ratio,
                              Array<Real>& dt_min, Array<Real>& dt_level,
                              Real stop_time, int post_regrid_flag);

    //
    // Print information about energy budget.
    //
    void do_energy_diagnostics();

    //
    // Do work after timestep().
    //
    virtual void post_timestep(int iteration);

    //
    // Contains operations to be done only after a full coarse timestep.
    //
    virtual void postCoarseTimeStep(Real cumtime);

    //
    // Do work after `regrid()`.
    //
    virtual void post_regrid(int lbase, int new_finest);

    //
    // Do work after a `restart()`.
    //
    virtual void post_restart();

    //
    // Do work after `init()`.
    //
    virtual void post_init(Real stop_time);

    //
    // Error estimation for regridding.
    //
    virtual void errorEst(TagBoxArray& tb, int clearval, int tagval, Real time,
                          int n_error_buf=0, int ngrow=0);

    //
    // Called in grid_places after other tagging routines to modify
    //   the list of tagged points
    //
    virtual void manual_tags_placement (TagBoxArray&    tags,
                                        Array<IntVect>& bf_lev);

    // Returns a MultiFab containing the derived data for this level. The user
    // is responsible for deleting this pointer when done with it. If
    // `ngrow` > 0 the MultiFab is built on the appropriately grown BoxArray.
    MultiFab* derive(const std::string& name, Real time, int ngrow);

    // This version of `derive()` fills the dcomp'th component of mf with the
    // derived quantity.
    void derive(const std::string& name, Real time, MultiFab& mf, int dcomp);

    static int Do_Hydro();
    static int num_grow();

    // Synchronize (rho e) and (rho E) so they are consistent with each other
    void reset_internal_energy(MultiFab& State, MultiFab& DiagEOS);

    void compute_new_temp();

    void compute_rho_temp(Real& rho_T_avg, Real& T_avg, Real& T_meanrho);

    void get_old_source(Real old_time, Real dt, MultiFab& Rhs);
    void get_new_source(Real old_time, Real new_time, Real dt, MultiFab& Rhs);

    Real vol_weight_sum(const std::string& name, Real time, bool masked);
    Real vol_weight_sum(MultiFab& mf, bool masked);

#ifdef DER_AMOM
    void vol_weight_sum_angmom(Real* angsum1, Real* angsum2, Real* angsum3, const Real time);
#endif
    Real vol_weight_squared_sum(const std::string& name, Real time);
    Real vol_weight_squared_sum_level(const std::string& name, Real time);

#ifdef AUX_UPDATE
    void advance_aux(Real time, Real dt);
#endif

    static Array<NyxParticleContainerBase*>& theActiveParticles();
    static Array<NyxParticleContainerBase*>& theVirtualParticles();
    static Array<NyxParticleContainerBase*>& theGhostParticles();

    static DarkMatterParticleContainer* theDMPC();
    static DarkMatterParticleContainer* theVirtPC();
    static DarkMatterParticleContainer* theGhostPC();

    static StellarParticleContainer* theSPC();
    static StellarParticleContainer* theVirtSPC();
    static StellarParticleContainer* theGhostSPC();

    static AGNParticleContainer* theAPC();
    static AGNParticleContainer* theVirtAPC();
    static AGNParticleContainer* theGhostAPC();

    static NeutrinoParticleContainer* theNPC();
    static NeutrinoParticleContainer* theVirtNPC();
    static NeutrinoParticleContainer* theGhostNPC();

    static int NUM_STATE;
    static int Density, Xmom, Ymom, Zmom, Eden, Eint;
    static int Temp_comp, Ne_comp;

    static int FirstSpec, FirstAux, FirstAdv;
    static int NumSpec, NumAux, NumAdv;

    static int strict_subcycling;

    static int init_with_sph_particles;

    static std::string particle_plotfile_format;

    //
    // This MultiFab is used for the level coarser than this level to mask out
    // this level.  We only build this when it is needed.
    // This MultiFab has to live on this level even though it is at the resolution
    // of the next coarser level because it must be updated whenever this level changes.
    //
    MultiFab* fine_mask;
    MultiFab* build_fine_mask();

    static int forceParticleRedist;  // ---- for dynamic sidecars
    static int nSidecarProcs;

    virtual void AddProcsToComp(Amr *aptr, int nSidecarProcs, int prevSidecarProcs,
                                int ioProcNumSCS, int ioProcNumAll, int scsMyId,
				MPI_Comm scsComm);
    virtual void NyxParticlesAddProcsToComp(Amr *aptr, int nSidecarProcs, int prevSidecarProcs,
                                            int ioProcNumSCS, int ioProcNumAll, int scsMyId,
				            MPI_Comm scsComm);
    static void InitErrorList();
    static void InitDeriveList();

protected:

    //
    // Number of grow cells required for particle motion.
    //
    int grav_n_grow;

    //
    // Initialize the network.
    //
    static void network_init();

    static void read_params();

    Nyx& get_level(int lev);

#ifndef NO_HYDRO
    FluxRegister& get_flux_reg();
    FluxRegister& get_flux_reg(int lev);

    void reflux();

    void enforce_nonnegative_species(MultiFab& S_new);
    void enforce_consistent_e(MultiFab& S);
#endif

    void average_down();
    void average_down(int state_indx);

    void build_metrics();

#ifndef NO_HYDRO
    virtual void sum_integrated_quantities();

    void compute_average_density();
    void compute_average_temperature(Real& average_temperature);
    void compute_average_species(int nspec, int naux, Array<Real>& average_species);

    void set_small_values();
#endif

    void write_info();

#ifndef NO_HYDRO
    FluxRegister* flux_reg;
#endif

    //
    // Static data members.
    //
    static bool dump_old;
    static int verbose;
    static int show_timings;
    static Real cfl;
    static Real init_shrink;
    static Real change_max;
    static int do_reflux;
    static ErrorList err_list;
    static BCRec phys_bc;
    static int NUM_GROW;

    static int nsteps_from_plotfile;

    static int allow_untagging;
    static int use_const_species;
    static int normalize_species;
    static int do_special_tagging;
    static int ppm_type;
    static int ppm_reference;
    static int ppm_flatten_before_integrals;
    static int use_colglaz;
    static int use_flattening;
    static int corner_coupling;
    static int version_2;

    static int use_exact_gravity;

    static Real small_dens;
    static Real small_temp;
    static Real gamma;

    static Real  h_species;
    static Real he_species;

    static bool do_dm_particles;

    // How do we want to initialize the particles?
    // Must be "Random", "Cosmological" or "AsciiFile"
    static std::string particle_init_type;

    // How do we want to move the particles?
    // Must be "Random" or "Gravitational"
    static std::string particle_move_type;

    // These control random initialization
    static bool particle_initrandom_serialize;
    static long particle_initrandom_count;
    static Real particle_initrandom_mass;
    static int particle_initrandom_iseed;

    static IntVect Nrep;  // how many times the initial conditions are replicated in each direction

    static Array<Real> plot_z_values;      // These are the value of "z" at which to dump plotfiles.
    static Array<Real> analysis_z_values;  // These are the value of "z" at which to perform analysis

    bool FillPatchedOldState_ok;

    // permits hydro to be turned on and off for running pure rad problems:
    static int do_hydro;

    // permits gravity calculation to be turned on and off
    static int do_grav;

    // if true, define an additional source term
    static int add_ext_src;

    // specifies the heating/cooling source term
    static int heat_cool_type;

    // if true , incorporate the source term through Strang-splitting
    // if false, incorporate the source term through predictor-corrector methodology
    static int strang_split;

#ifdef GRAVITY
    // There can be only one Gravity object, it covers all levels:
    static class Gravity *gravity;
#endif

    // Average value of gas, dark matter, neutrino and total density, computed every coarse timestep
    static Real average_gas_density;
    static Real average_dm_density;
    static Real average_neutr_density;
    static Real average_total_density;


  // for keeping track of the amount of CPU time used -- this will persist
  // after restarts
  static Real      previousCPUTimeUsed;
  static Real      startCPUTime;

  static Real getCPUTime();


};

// time step interval for finding halos
extern int reeber_int;
// time step interval for doing Gimlet post-processing
extern int gimlet_int;

//
// Inlines.
//
inline
int
Nyx::Do_Hydro()
{
    return do_hydro;
}

inline
int
Nyx::num_grow()
{
    return NUM_GROW;
}

inline
Nyx&
Nyx::get_level(int level)
{
    return *(Nyx *) &parent->getLevel(level);
}

#ifndef NO_HYDRO
inline
FluxRegister&
Nyx::get_flux_reg()
{
    BL_ASSERT(flux_reg);
    return *flux_reg;
}

inline
FluxRegister&
Nyx::get_flux_reg(int level)
{
    return get_level(level).get_flux_reg();
}
#endif // NO_HYDRO

#endif /*_Nyx_H_*/
