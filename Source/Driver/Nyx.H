
#ifndef _Nyx_H_
#define _Nyx_H_

#include <AMReX_BC_TYPES.H>
#include <AMReX_AmrLevel.H>
#include <AMReX_ErrorList.H>
#include <AMReX_FluxRegister.H>

#ifdef AMREX_PARTICLES
#include <NyxParticleContainer.H>
#include <DarkMatterParticleContainer.H>
#ifdef AGN
#include <AGNParticleContainer.H>
#endif
#ifdef NEUTRINO_PARTICLES
#include <NeutrinoParticleContainer.H>
#endif
#else
#include <AMReX_ParmParse.H>
#endif

#ifdef REEBER
#ifndef REEBER_HIST
#include <reeber.H>
#endif
#endif

#ifdef AMREX_USE_CONDUIT
#include <AMReX_Conduit_Blueprint.H>
using namespace amrex;
using namespace conduit;

#include <ascent.hpp>
using namespace ascent;
#endif

#include <iostream>

#ifdef BL_HDF5
#include <hdf5.h>
#endif

using std::istream;
using std::ostream;
#ifdef AMREX_PARTICLES
typedef NyxParticleContainer<1+AMREX_SPACEDIM+3> StellarParticleContainer;
#endif
#ifdef NO_HYDRO
enum StateType {
    PhiGrav_Type = 0,
    Gravity_Type,
    NUM_STATE_TYPE
};
static constexpr int max_prob_param = 5;
#endif

#ifndef NO_HYDRO
enum StateType {
    State_Type = 0,
    DiagEOS_Type,
    PhiGrav_Type,
    Gravity_Type,
#ifdef SDC
    SDC_IR_Type,
#endif
    NUM_STATE_TYPE
};
static constexpr int max_prob_param = 20;
#endif

//
// AmrLevel-derived class for hyperbolic conservation equations for stellar
// media
//

class Nyx
    :
    public amrex::AmrLevel
{
public:
    //
    //Default constructor.  Builds invalid object.
    //
    Nyx();

    //
    //The basic constructor.
    //
    Nyx(amrex::Amr& papa, int lev, const amrex::Geometry& level_geom,
        const amrex::BoxArray& bl, const amrex::DistributionMapping& dm,
        amrex::Real time);

    //
    //The destructor.
    //
    virtual ~Nyx() override;

    //
    //Restart from a checkpoint file.
    //
    virtual void restart(amrex::Amr& papa, istream& is, bool b_read_special=false) override;

    //
    //Call amrex::AmrLevel::checkPoint and then add radiation info
    //
    virtual void checkPoint(const std::string& dir, std::ostream& os,
                            amrex::VisMF::How how, bool dump_old) override;
    virtual void checkPointPre(const std::string& dir, std::ostream& os) override;
    virtual void checkPointPost(const std::string& dir, std::ostream& os) override;

    // A string written as the first item in `write_plot_file()` at level zero.
    // It is so we can distinguish between different types of plot files. For
    // Nyx it has the form: Nyx-Vnnn.
    virtual std::string thePlotFileType() const override;

    virtual void setPlotVariables() override;

    //
    //Write a plotfile to specified directory.
    //
    virtual void writePlotFile(const std::string& dir, ostream& os, amrex::VisMF::How how) override;
#ifdef BL_HDF5
    virtual void writePlotFileHDF5(const std::string& dir, ostream& os, amrex::VisMF::How how) override;
#endif
    virtual void writePlotFilePre(const std::string& dir, ostream& os) override;
    virtual void writePlotFilePost(const std::string& dir, ostream& os) override;

    static void writeBuildInfo ();
    void writeJobInfo (const std::string& dir);

#ifdef AMREX_USE_CONDUIT
    void blueprint_check_point ();
#endif

    //
    //Write amrex::MultiFab as plot file
    //
    void writeMultiFabAsPlotFile(const std::string& pltfile, const amrex::MultiFab& mf, std::string componentName);
    //
    //Write all parameters into specified directory.
    //
    static int write_parameters_in_plotfile;
    virtual void write_parameter_file(const std::string& dir);
    static int write_skip_prepost;
    static int write_hdf5;

    static int runlog_precision;
    static int runlog_precision_terse;

    //
    //Define data descriptors.
    //
    static void variable_setup();
    static void hydro_setup();
    static void heatcool_setup();
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
    virtual void initData() override;

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
    // Initial redshift
    //
    static amrex::Real initial_z;

    //
    // Initial time in code units
    //
    static amrex::Real initial_time;

    //
    // End time in code units
    //
    static amrex::Real final_time;

    //
    // Final a -- used as stopping criterion if positive
    //
    static amrex::Real final_a;

    //
    // Final z -- used as stopping criterion if positive
    //
    static amrex::Real final_z;

    //
    // Relative change in a allowed in one timestep
    //
    static amrex::Real relative_max_change_a;

    //
    // Absolute change in a allowed in one timestep for fixed delta_a
    //
    static amrex::Real absolute_max_change_a;

    //
    // Positive number means use powers of 2 binning for relative dt
    //
    static amrex::Real dt_binpow;

    //
    // Old and new times at which "old_a" and "new_a" are defined.
    //
    static amrex::Real old_a_time;
    static amrex::Real new_a_time;

    //
    // "a" at old_a_time and new_a_time
    //
    static amrex::Real old_a;
    static amrex::Real new_a;

    //
    // comoving parameters
    //
    static amrex::Real comoving_OmB;
    static amrex::Real comoving_OmM;
    static amrex::Real comoving_OmR;
    static amrex::Real comoving_h;
    static         int comoving_type;

    //
    // Get the comoving coordinate "a"
    //
    amrex::Real get_comoving_a(amrex::Real time);
    static void integrate_time_given_a(const amrex::Real a0, const amrex::Real a1, amrex::Real& dt);
    void integrate_comoving_a(const amrex::Real old_a, amrex::Real& new_a, const amrex::Real dt);
    void integrate_comoving_a(amrex::Real time, amrex::Real dt);

    //
    // Time step control based on "a" not growing too fast
    //
    void integrate_comoving_a_to_a(const amrex::Real old_a, const amrex::Real a_value, amrex::Real& dt);
    void integrate_comoving_a_to_z(const amrex::Real old_a, const amrex::Real z_value, amrex::Real& dt);
    void est_maxdt_comoving_a(const amrex::Real old_a, amrex::Real & dt);
    void est_lindt_comoving_a(const amrex::Real old_a, const amrex::Real new_a, amrex::Real& dt);
    void estdt_comoving_a(const amrex::Real old_a,amrex::Real& new_a,amrex::Real& dt,const amrex::Real change_allowed,const amrex::Real fixed_da,const amrex::Real final_a,int& dt_modified);
    void enforce_percent_change(const amrex::Real old_a,amrex::Real& new_a,amrex::Real& dt,const amrex::Real change_allowed);
    void enforce_final_a(const amrex::Real old_a,amrex::Real& new_a,amrex::Real& dt,const amrex::Real final_a);
    void comoving_est_time_step(amrex::Real& cur_time, amrex::Real& est_dt);

    //
    // Time step control based on "z" not passing one of the specified plot_z_values
    //
    void plot_z_est_time_step(amrex::Real& est_dt, bool& dt_changed);

    //
    // Time step control based on "z" not passing one of the specified analysis_z_values
    //
    void analysis_z_est_time_step(amrex::Real& est_dt, bool& dt_changed);

    //
    // How to initialize "a" at restart (from checkpoint or plotfile)
    //
    void comoving_a_post_restart(const std::string& restart_file);

    //
    // How to initialize typical values at restart (from checkpoint or plotfile)
    //
    void typical_values_post_restart(const std::string& restart_file);

    //
    // Initialize the zhi component of the EOS from a binary input file
    //
    void init_zhi ();

    //
    // Initialize the Santa Barbara problem
    //
    void init_santa_barbara(int init_sb_vels);

    void check_initial_species();
    void init_e_from_T(const amrex::Real& a);

    //
    // Initialize from MUSIC
    //
    void initcosmo();
    void icReadAndPrepareFab(std::string mfDirName, int nghost, amrex::MultiFab &mf);

    //
    // Write forcing spectrum in checkpoint directories
    //
    void forcing_check_point (const std::string& dir);

    //
    // Restore spectrum at restart
    //
    void forcing_post_restart(const std::string& restart_file);

    void integrate_state_force(
            amrex::Box const& bx,
            amrex::Array4<amrex::Real> const& state,
            amrex::Array4<amrex::Real> const& diag_eos,
            const amrex::Real* dx,
            const amrex::Real time,
            const amrex::Real a,
            const amrex::Real half_dt);

        void ext_src_force(
            amrex::Box const& bx,
            amrex::Array4<const amrex::Real> const& old_state,
            amrex::Array4<const amrex::Real> const& new_state,
            amrex::Array4<const amrex::Real> const& old_diag,
            amrex::Array4<amrex::Real> const& new_diag,
            amrex::Array4<amrex::Real> const& src,
            const amrex::Real* problo,
            const amrex::Real* dx,
            const amrex::Real time,
            const amrex::Real z,
            const amrex::Real dt);

    //
    // Read particle-related inputs
    //
    static void read_particle_params();

    //
    // Read initialization-related inputs
    //
    static void read_init_params();

    //
    // Read hydro-related inputs
    //
#ifndef NO_HYDRO
    static void read_hydro_params();
#endif

    //
    // Write particles in checkpoint directories
    //
    void particle_check_point(const std::string& dir);

#ifdef AMREX_PARTICLES
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
    // override for sanity
    void particle_redistribute(int lbase = 0, bool init = false) override;

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
    void setup_ghost_particles(int ngrow);

    //
    // Remove ghost particles (for this level) if necessary
    //
    void remove_ghost_particles();

    //
    // Time step control based on particles
    //
    void particle_est_time_step(amrex::Real& est_dt);

    //
    // Default verbosity of Particle class
    //
    static int particle_verbose;

    //
    // Default cfl of particles in Particle class
    //
    static amrex::Real particle_cfl;
#ifdef NEUTRINO_PARTICLES
    static amrex::Real neutrino_cfl;
#endif

    //
    // Shall we write the initial single-level particle density into a multifab
    //   called "ParticleDensity"?
    //
    static int write_particle_density_at_init;
#endif

    //
    // Derived quantities associated with particles
    //
    std::unique_ptr<amrex::MultiFab> particle_derive (const std::string& name, amrex::Real time, int ngrow);

    //
    //Set time levels of state data.
    //
    virtual void setTimeLevel (amrex::Real time, amrex::Real dt_old, amrex::Real dt_new) override;

    //
    //Initialize data on this level from another Nyx (during regrid).
    //
    virtual void init(amrex::AmrLevel& old) override;

    //
    // Initialize data on this level after regridding if old level did not
    // previously exist
    //
    virtual void init() override;

    //
    // Proceed with next timestep?
    //
    virtual int okToContinue() override;

    //
    // Tell amrex::Amr to write a plotfile now
    //
    bool writePlotNow() override;

    //
    // Tell Nyx to do analysis now
    //
    bool doAnalysisNow();

    int updateInSitu();

    //
    // Advance grids at this level in time.
    //
    virtual amrex::Real advance(amrex::Real time, amrex::Real dt, int iteration, int ncycle) override;

    amrex::Real advance_hydro(amrex::Real time, amrex::Real dt, int iteration, int ncycle);
    amrex::Real advance_heatcool(amrex::Real time, amrex::Real dt, int iteration, int ncycle);
    amrex::Real advance_no_hydro(amrex::Real time, amrex::Real dt, int iteration, int ncycle);
    amrex::Real advance_hydro_plus_particles(amrex::Real time, amrex::Real dt, int iteration,
                                      int ncycle);

    void strang_hydro(amrex::Real time, amrex::Real dt, amrex::Real a_old, amrex::Real a_new);
#ifdef SDC
    void    sdc_hydro(amrex::Real time, amrex::Real dt, amrex::Real a_old, amrex::Real a_new);
#endif

    void correct_gsrc(int lev, amrex::Real time, amrex::Real prev_time, amrex::Real cur_time,
                               amrex::Real dt);

///
/// this constructs the hydrodynamic source (essentially the flux
/// divergence) using the CTU framework for unsplit hydrodynamics
///
/// @param time     current time
/// @param dt       timestep
///
    void construct_ctu_hydro_source(amrex::Real time, amrex::Real dt, amrex::Real a_old, amrex::Real a_new,
                               amrex::MultiFab& S_border, amrex::MultiFab& D_border,
                               amrex::MultiFab& ext_src_old, amrex::MultiFab& hydro_src,
                               amrex::MultiFab& grav,
                               bool init_flux_register, bool add_to_flux_register);

///
/// Hydrodynamic (and radiation) fluxes.
///
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > fluxes;
///
/// Return the n'th fluxes MultiFab.
///
/// @param dir  direction in which to find the fluxes along
///
    amrex::MultiFab& Fluxes (int dir);

    amrex::Vector<std::unique_ptr<amrex::MultiFab> > mass_fluxes;

    void construct_hydro_source(const amrex::MultiFab& S,
                                      amrex::MultiFab& sources_for_hydro,
                                      amrex::MultiFab& hydro_source,
                                      amrex::MultiFab& grav,
                                      amrex::Real a_old,
                                      amrex::Real a_new,
                                      amrex::Real dt,
                                      bool init_flux_register, bool add_to_flux_register);

    void update_state_with_sources( amrex::MultiFab& S_old      , amrex::MultiFab& S_new,
                                    amrex::MultiFab& ext_src_old, amrex::MultiFab& hydro_src,
                                    amrex::MultiFab& grav       , 
#ifdef SDC
                                    amrex::MultiFab& reset_e_src,
#endif
                                    amrex::Real dt, amrex::Real a_old, amrex::Real a_new);

    void enforce_minimum_density  ( amrex::MultiFab& S_old      , amrex::MultiFab& S_new, 
#ifdef SDC
                                    amrex::MultiFab& hydro_source,
                                    amrex::MultiFab& reset_e_src,
#endif
                                    amrex::Real a_new);
#ifdef SDC
    void enforce_minimum_density_cons  ( amrex::MultiFab& S_old      , amrex::MultiFab& S_new, amrex::MultiFab& reset_e_src);
#else
    void enforce_minimum_density_cons  ( amrex::MultiFab& S_old      , amrex::MultiFab& S_new);
#endif
    void enforce_minimum_density_floor ( amrex::MultiFab& S_new, amrex::Real a_new);

    void strang_first_step  (amrex::Real time, amrex::Real dt,  amrex::MultiFab& state, amrex::MultiFab&  dstate);
    void strang_second_step (amrex::Real time, amrex::Real dt,  amrex::MultiFab& state, amrex::MultiFab&  dstate);

    void sdc_reactions ( amrex::MultiFab& state_old, amrex::MultiFab& state_new, amrex::MultiFab&  dstate,
                         amrex::MultiFab& hydro_src, amrex::MultiFab& IR, amrex::MultiFab& reset_e_src,
                         amrex::Real dt, amrex::Real a_old, amrex::Real a_new,
                         int sdc_iter);

  int integrate_state_vec(amrex::MultiFab &state,   amrex::MultiFab &diag_eos, const amrex::Real& a, const amrex::Real& delta_time);
   int integrate_state_grownvec(amrex::MultiFab &state,   amrex::MultiFab &diag_eos, const amrex::Real& a, const amrex::Real& delta_time);
  int integrate_state_vec_mfin(amrex::Array4<amrex::Real>const& state4,   amrex::Array4<amrex::Real>const& diag_eos4,const  amrex::Box& tbx,  const amrex::Real& a, const amrex::Real& delta_time, long int& old_max_steps, long int& new_max_steps);

int integrate_state_struct
  (amrex::MultiFab &S_old,
   amrex::MultiFab &S_new,
   amrex::MultiFab &D_old,
   amrex::MultiFab &hydro_src,
   amrex::MultiFab &IR,
   amrex::MultiFab &reset_src,
   const amrex::Real& a, const amrex::Real& a_end,
   const amrex::Real& delta_time,
   const int sdc_iter);

  int integrate_state_struct_mfin
  (amrex::Array4<amrex::Real> const& state4,
   amrex::Array4<amrex::Real> const& diag_eos4,
   amrex::Array4<amrex::Real> const& state_n4,
   amrex::Array4<amrex::Real> const& hydro_src4,
   amrex::Array4<amrex::Real> const& reset_src4,
   amrex::Array4<amrex::Real> const& IR4,
#ifdef SAVE_REACT
   amrex::Array4<amrex::Real> const& react_in_arr,
   amrex::Array4<amrex::Real> const& react_out_arr,
   amrex::Array4<amrex::Real> const& react_out_work_arr,
#endif
   const amrex::Box& tbx,
   const amrex::Real& a, const amrex::Real& a_end,
   const amrex::Real& delta_time,
   long int& old_max_steps, long int& new_max_steps,
   const int sdc_iter);

  int integrate_state_cell(amrex::MultiFab &state,   amrex::MultiFab &diag_eos, const amrex::Real& a, const amrex::Real& delta_time);
   int integrate_state_growncell(amrex::MultiFab &state,   amrex::MultiFab &diag_eos, const amrex::Real& a, const amrex::Real& delta_time);

    amrex::Real advance_particles_only (amrex::Real time, amrex::Real dt, int iteration, int ncycle);

    void moveKickDriftExact(amrex::Real dt);
    void moveKickExact(amrex::Real dt);
    void time_center_source_terms(amrex::MultiFab& S_new, amrex::MultiFab& ext_src_old,
                                  amrex::MultiFab& ext_src_new, amrex::Real dt);

    void conserved_to_primitive(amrex::MultiFab& state);
    void primitive_to_conserved(amrex::MultiFab& state);

#ifdef REEBER
    void halo_find(amrex::Real dt,        amrex::Vector<Halo>& reeber_halos);
    void halo_print(amrex::Vector<Halo>& reeber_halos);
#endif
    void agn_halo_find(amrex::Real dt);
    void agn_halo_merge();
    void agn_halo_accrete(amrex::Real dt);

#ifdef REEBER
    void runReeberAnalysis(amrex::Vector<amrex::MultiFab*>& new_state,
                                       amrex::Vector<std::unique_ptr<amrex::MultiFab> >& particle_mf,
                           const amrex::Geometry Geom, const amrex::Vector<amrex::IntVect>& level_refinements,
                           int nStep, bool do_analysis, std::vector<Halo>& reeber_halos);
#endif

    void Lya_statistics();

    //
    // Estimate time step.
    //
    amrex::Real est_time_step(amrex::Real dt_old);

    //
    // Compute initial time step.
    //
    amrex::Real initial_time_step();

    //
    // Compute initial `dt'.
    //
    virtual void computeInitialDt(int finest_level, int sub_cycle,
                                  amrex::Vector<int>& n_cycle,
                                  const amrex::Vector<amrex::IntVect>& ref_ratio,
                                  amrex::Vector<amrex::Real>& dt_level, amrex::Real stop_time) override;
    //
    // Compute new `dt'.
    //
    virtual void computeNewDt(int finest_level, int sub_cycle,
                              amrex::Vector<int>& n_cycle,
                              const amrex::Vector<amrex::IntVect>& ref_ratio,
                              amrex::Vector<amrex::Real>& dt_min, amrex::Vector<amrex::Real>& dt_level,
                              amrex::Real stop_time, int post_regrid_flag) override;

    //
    // Print information about energy budget.
    //
    void do_energy_diagnostics();

    //
    // Do work after timestep().
    //
    virtual void post_timestep(int iteration) override;

    //
    // Contains operations to be done only after a full coarse timestep.
    //
    virtual void postCoarseTimeStep(amrex::Real cumtime) override;

    //
    // Do work after `regrid()`.
    //
    virtual void post_regrid(int lbase, int new_finest) override;

    //
    // Do work after a `restart()`.
    //
    virtual void post_restart() override;

    //
    // Do work after `init()`.
    //
    virtual void post_init(amrex::Real stop_time) override;

    //
    // Error estimation for regridding.
    //
    virtual void errorEst(amrex::TagBoxArray& tb, int clearval, int tagval, amrex::Real time,
                          int n_error_buf=0, int ngrow=0) override;

    //
    // Error estimation for regridding.
    //
    virtual int WorkEstType () override { return 0; }

    //
    // Called in grid_places after other tagging routines to modify
    //   the list of tagged points
    //
    virtual void manual_tags_placement (amrex::TagBoxArray&    tags,
                                        const amrex::Vector<amrex::IntVect>& bf_lev) override;

    // Returns a amrex::MultiFab containing the derived data for this level. The user
    // is responsible for deleting this pointer when done with it. If
    // `ngrow` > 0 the amrex::MultiFab is built on the appropriately grown amrex::BoxArray.
    std::unique_ptr<amrex::MultiFab> derive(const std::string& name, amrex::Real time, int ngrow) override;

    // This version of `derive()` fills the dcomp'th component of mf with the
    // derived quantity.
    void derive(const std::string& name, amrex::Real time, amrex::MultiFab& mf, int dcomp) override;

    static int Do_Hydro();

#ifndef CONST_SPECIES
    void enforce_nonnegative_species(amrex::MultiFab& S_new);
#endif
    void enforce_consistent_e(amrex::MultiFab& S);

    // Synchronize (rho e) and (rho E) so they are consistent with each other
    void reset_internal_energy(amrex::MultiFab& State, amrex::MultiFab& DiagEOS, amrex::MultiFab& reset_e_src);

    // Synchronize (rho e) and (rho E) so they are consistent with each other without storing the change
    void reset_internal_energy_nostore(amrex::MultiFab& State, amrex::MultiFab& DiagEOS);

    void reset_internal_energy_interp(amrex::MultiFab& State, amrex::MultiFab& DiagEOS, amrex::MultiFab& reset_e_src);

    // Note: this no longer includes the call to reset_internal_energy
    void compute_new_temp(amrex::MultiFab& S_new, amrex::MultiFab& D_new);

    void compute_rho_temp(amrex::Real& rho_T_avg, amrex::Real& T_avg, amrex::Real& Tinv_avg, amrex::Real& T_meanrho);
    void compute_gas_fractions(amrex::Real T_cut, amrex::Real rho_cut,
                               amrex::Real& whim_mass_frac, amrex::Real& whim_vol_frac,
                               amrex::Real& hh_mass_frac, amrex::Real& hh_vol_frac,
                               amrex::Real& igm_mass_frac, amrex::Real& igm_vol_frac);

    void get_old_source(amrex::Real old_time, amrex::Real dt, amrex::MultiFab& Rhs);
    void get_new_source(amrex::Real old_time, amrex::Real new_time, amrex::Real dt, amrex::MultiFab& Rhs);

    amrex::Real vol_weight_sum              (const std::string& name, amrex::Real time, bool masked);
    amrex::Real vol_weight_sum              (amrex::MultiFab& mf, bool masked);
    amrex::Real vol_weight_squared_sum      (const std::string& name, amrex::Real time);
    amrex::Real vol_weight_squared_sum_level(const std::string& name, amrex::Real time);

#ifdef AUX_UPDATE
    void advance_aux(amrex::Real time, amrex::Real dt);
#endif
#ifdef AMREX_PARTICLES
    static amrex::Vector<NyxParticleContainerBase*>& theActiveParticles();
    static amrex::Vector<NyxParticleContainerBase*>& theVirtualParticles();
    static amrex::Vector<NyxParticleContainerBase*>& theGhostParticles();

    static DarkMatterParticleContainer* theDMPC();
    static DarkMatterParticleContainer* theVirtPC();
    static DarkMatterParticleContainer* theGhostPC();

    static StellarParticleContainer* theSPC();
    static StellarParticleContainer* theVirtSPC();
    static StellarParticleContainer* theGhostSPC();

#ifdef AGN
    static AGNParticleContainer* theAPC();
    static AGNParticleContainer* theVirtAPC();
    static AGNParticleContainer* theGhostAPC();
#endif

#ifdef NEUTRINO_PARTICLES
  #ifdef NEUTRINO_DARK_PARTICLES
    static DarkMatterParticleContainer* theNPC();
    static DarkMatterParticleContainer* theVirtNPC();
    static DarkMatterParticleContainer* theGhostNPC();
  #else
    static NeutrinoParticleContainer* theNPC();
    static NeutrinoParticleContainer* theVirtNPC();
    static NeutrinoParticleContainer* theGhostNPC();
  #endif
#endif
#endif

    static int NUM_STATE;

    static int Zhi_comp;

    static int NumSpec;
    static amrex::Real gamma;
    static amrex::Real h_species;
    static amrex::Real he_species;

    static int init_with_sph_particles;

    static std::string particle_plotfile_format;

    // Average value of gas, dark matter, neutrino and total density, computed every coarse timestep
    static amrex::Real average_gas_density;
    static amrex::Real average_dm_density;
    static amrex::Real average_neutr_density;
    static amrex::Real average_total_density;

    static amrex::Real tagging_base;

    // specifies the memory priority
    static int reuse_mlpoisson;
  
    //
    // This amrex::MultiFab is used for the level coarser than this level to mask out
    // this level.  We only build this when it is needed.
    // This amrex::MultiFab has to live on this level even though it is at the resolution
    // of the next coarser level because it must be updated whenever this level changes.
    //
    amrex::MultiFab* fine_mask;
    amrex::MultiFab* build_fine_mask();

    static void InitErrorList();
    static void InitDeriveList();

    //! Get the level directory names
    void LevelDirectoryNames (const std::string &dir,
                              const std::string &secondDir,  // ---- probably DM or AGN
                              std::string &LevelDir,
                              std::string &FullPath);
    //! Create the Level_ and other directories in checkpoint and plot files
    virtual void CreateLevelDirectory (const std::string &dir) override;

protected:

    static void read_params();

    Nyx& get_level(int lev);

    std::string retrieveDM();
#ifdef AGN
    std::string retrieveAGN();
#endif
#ifdef NEUTRINO_DARK_PARTICLES
    std::string retrieveNPC();
#endif

#ifndef NO_HYDRO
    amrex::FluxRegister& get_flux_reg();
    amrex::FluxRegister& get_flux_reg(int lev);

    void reflux();

#endif

    void average_down();
    void average_down(int state_indx);

    void build_metrics();

#ifndef NO_HYDRO
    virtual void sum_integrated_quantities();

    void compute_average_density();
    void compute_average_temperature(amrex::Real& average_temperature);
    void compute_average_species(amrex::Vector<amrex::Real>& average_species);

    void set_small_values();

#endif

    void set_small_values_given_average (amrex::Real average_dens, amrex::Real average_temp, amrex::Real a,
                                         amrex::Real & small_dens_inout, amrex::Real & small_temp_inout,
                                         amrex::Real &small_pres_inout, amrex::Real gamma_minus_1,
                                         amrex::Real h_species);

    void write_info();

#ifndef NO_HYDRO
    amrex::FluxRegister* flux_reg;
#endif

    //
    // Static data members.
    //
    static bool dump_old;
    static int verbose;
    static amrex::Real cfl;
    static amrex::Real init_shrink;
    static amrex::Real change_max;
    static int do_reflux;
    static amrex::ErrorList err_list;
    static amrex::Vector<amrex::AMRErrorTag> errtags;
    static amrex::BCRec phys_bc;

    static int nsteps_from_plotfile;

    static int nghost_state;
    static int ppm_type;
    static int strang_grown_box;

    static std::string enforce_min_density_type;

    static int num_particle_ghosts;

    static amrex::Real small_dens;
    static amrex::Real small_pres;
    static amrex::Real small_temp;
    static amrex::Real small;

    static amrex::Real large_temp;

    static bool do_dm_particles;

    // How do we want to initialize the particles?
    // Must be "Random", "Cosmological" or "AsciiFile"
    static std::string particle_init_type;

    // These control random initialization
    static bool particle_initrandom_serialize;
    static long particle_initrandom_count;
    static long particle_initrandom_count_per_box;
    static amrex::Real particle_initrandom_mass;
    static amrex::Real particle_initrandom_mass_total;
    static int particle_initrandom_iseed;
    static int particle_skip_factor;
    static amrex::Real particle_inituniform_mass;
    static amrex::Real particle_inituniform_vx;
    static amrex::Real particle_inituniform_vy;
    static amrex::Real particle_inituniform_vz;

    // This controls the LaunchSafeGuard for particle initialization
    static int particle_launch_ics;

    static amrex::Vector<amrex::Real> plot_z_values;      // These are the value of "z" at which to dump plotfiles.
    static amrex::Vector<amrex::Real> analysis_z_values;  // These are the value of "z" at which to perform analysis
    static int insitu_int;
    static int insitu_start;

    static int load_balance_int;
    static amrex::Real load_balance_start_z;
    static int load_balance_wgt_strategy;
    static int load_balance_wgt_nmax;
    static int load_balance_strategy;

    bool FillPatchedOldState_ok;

    // permits hydro to be turned on and off for running pure rad problems:
    static int do_hydro;

    // permits gravity calculation to be turned on and off
    static int do_grav;

    // if true, define an additional source term
    static int add_ext_src;

    // specifies the heating/cooling source term
    static int heat_cool_type;

    // specifies the reduction type for snudials vectors
    static int sundials_atomic_reductions;

    // specifies the allocation type for sundials vector
    static int sundials_alloc_type;

    // specifies the allocation type for sundials vector
    static int use_typical_steps;

    // if true, use non-negative constraint on energy
    static int use_sundials_constraint;

    // if true, use non-negative constraint on energy
    static int use_sundials_fused;

    // specifies whether we allow tiling when calling the SUNDIALS integrator
    static bool sundials_use_tiling;

    // specifies tolerances when calling the SUNDIALS integrator
    static amrex::Real sundials_reltol;
    static amrex::Real sundials_abstol;

    // specifies tile_size of advance_hydro for unsplit Godonuv hydro advance
    static amrex::IntVect     hydro_tile_size;

    // specifies tile_size of integrate_state for heating-cooling reactions
    static amrex::IntVect      sundials_tile_size;

    // specifies the memory priority
    static int minimize_memory;

    // specifies the ShrinkToFit priority
    static int shrink_to_fit;

    // specifies inhomogeneous reionization type
    static int inhomo_reion;
    static std::string inhomo_zhi_file;
    static int inhomo_grid;

    // permits forcing to be switched on and off
    static int do_forcing;

    // if true , incorporate the source term through Strang-splitting
    static int strang_split;

#ifdef SDC
    // if true , use SDC to couple hydro and reactions
    static int sdc_split;
    static int strang_restart_from_sdc;
#endif

    // There can be only one Gravity object, it covers all levels:
    static class Gravity *gravity;

    // There can be only one forcing object (Fourier space):
    static class StochasticForcing *forcing;

#ifdef REEBER
  //
  // Threshold for halo to create SMBH
  //
  static amrex::Real mass_halo_min;

  //
  // Seed mass of SMBH
  //
  static amrex::Real mass_seed;
#endif

  // Previous maximum number of steps for sundials
  static long int old_max_sundials_steps;
  static long int new_max_sundials_steps;

  // for keeping track of the amount of CPU time used -- this will persist
  // after restarts
  static amrex::Real      previousCPUTimeUsed;
  static amrex::Real      startCPUTime;

  static amrex::Real getCPUTime();

};

#ifndef NO_HYDRO
#include <reset_internal_e.H>
#endif

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
Nyx&
Nyx::get_level(int my_level)
{
    return *(Nyx *) &parent->getLevel(my_level);
}

#ifndef NO_HYDRO
inline
amrex::FluxRegister&
Nyx::get_flux_reg()
{
    BL_ASSERT(flux_reg);
    return *flux_reg;
}

inline
amrex::FluxRegister&
Nyx::get_flux_reg(int my_level)
{
    return get_level(my_level).get_flux_reg();
}
#endif // NO_HYDRO

#endif /*_Nyx_H_*/
