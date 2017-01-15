#ifndef _Gravity_H_
#define _Gravity_H_

#include <AMReX_AmrLevel.H>
#include <AMReX_MacBndry.H>
#include <AMReX_FluxRegister.H>
#include <AMReX_Particles.H>

class Gravity {

public:

    Gravity(amrex::Amr* Parent, int _finest_level, amrex::BCRec* _phys_bc, int _Density);
    virtual ~Gravity();

    void read_params();

    void install_level(int level, amrex::AmrLevel* level_data_to_install);

    std::string get_gravity_type();
    Real get_const_grav ();
    int get_no_sync();
    int get_no_composite();

    void set_mass_offset(Real time);

    amrex::Array<amrex::MultiFab*> get_grad_phi_prev(int level);
    amrex::Array<amrex::MultiFab*> get_grad_phi_curr(int level);

    void plus_phi_curr(int level, amrex::MultiFab& addend);
    void plus_grad_phi_curr(int level, const amrex::Array<amrex::MultiFab*>& addend);

    // 
    // Moves curr phi data to prev and resets curr to 1.e50
    //
    void swap_time_levels(int level);

    // 
    // Use the F90 solver
    //
    Real solve_with_fmg (int crse_level, int fine_level,
                         const amrex::Array<amrex::MultiFab*> & phi,
                         const amrex::Array<amrex::MultiFab*> & rhs,
                         const amrex::Array<amrex::Array<amrex::MultiFab*> >& grad_phi,
                         const amrex::Array<amrex::MultiFab*>& res,
                         Real time);

    // 
    // Use the C++ solver in C_CellMG instead of the F90 solvers in F_MG.
    //
    void solve_with_Cpp(int level, amrex::MultiFab& phi, const amrex::Array<amrex::MultiFab*>& grad_phi,
                        amrex::MultiFab& rhs, Real tol, Real abs_tol);
    // 
    // Use the HPGMG solver
    //
#ifdef USEHPGMG
    void solve_with_HPGMG(int level, amrex::MultiFab& phi, const amrex::Array<amrex::MultiFab*>& grad_phi,
                        amrex::MultiFab& rhs, Real tol, Real abs_tol);
#endif
    void set_boundary  (amrex::BndryData& bd, amrex::MultiFab&  rhs, const Real* dx);

    void solve_for_old_phi(int level, amrex::MultiFab& phi, const amrex::Array<amrex::MultiFab*>& grad_phi,
                           int fill_interior);
    void solve_for_new_phi(int level, amrex::MultiFab& phi, const amrex::Array<amrex::MultiFab*>& grad_phi,
                           int fill_interior, int grav_n_grow = 1);
    void solve_for_phi(int level, amrex::MultiFab& Rhs, amrex::MultiFab& phi,
                       const amrex::Array<amrex::MultiFab*>& grad_phi, Real time,
                       int fill_interior);
    void solve_for_delta_phi(int crse_level, int fine_level, amrex::MultiFab& CrseRhs,
                             const amrex::Array<amrex::MultiFab*>& delta_phi,
                             const amrex::Array<amrex::Array<amrex::MultiFab*> >& grad_delta_phi);

    void gravity_sync(int crse_level, int fine_level, int iteration, int ncycle, 
                      const amrex::MultiFab& drho_and_drhoU, const amrex::MultiFab& dphi,
                      const amrex::Array<amrex::MultiFab*>& grad_delta_phi_cc);

    void multilevel_solve_for_old_phi(int level, int finest_level,
                                      int use_previous_phi_as_guess=0);
    void multilevel_solve_for_new_phi(int level, int finest_level,
                                      int use_previous_phi_as_guess=0);
    void multilevel_solve_for_phi(int level, int finest_level,
                                  int use_previous_phi_as_guess=0);
    void actual_multilevel_solve(int level, int finest_level,
                                 const amrex::Array<amrex::Array<amrex::MultiFab*> >& grad_phi, int is_new,
                                 int use_previous_phi_as_guess=0);

    void get_crse_grad_phi(int level, amrex::Array<std::unique_ptr<amrex::MultiFab> >& grad_phi_crse,
                           Real time);
    void get_crse_phi(int level, amrex::MultiFab& phi_crse, Real time);

    //
    // Sets phi_flux registers to 0
    //
    void zero_phi_flux_reg(int level);

    void get_old_grav_vector(int level, amrex::MultiFab& grav_vector, Real time);
    void get_new_grav_vector(int level, amrex::MultiFab& grav_vector, Real time);

    void average_fine_ec_onto_crse_ec(int level, int is_new);

    void add_to_fluxes(int level, int iteration, int ncycle);

    void reflux_phi(int level, amrex::MultiFab& dphi);

    void make_mg_bc();

    void set_dirichlet_bcs(int level, amrex::MultiFab* phi);

#ifdef CGRAV
    void make_prescribed_grav(int level, Real time, amrex::MultiFab& grav, int addToExisting);
#endif

    // Routine to duplicate Gravity class data onto sidecars
    virtual void AddProcsToComp(amrex::Amr *aptr, int level, amrex::AmrLevel *level_data_to_install,
                                int ioProcNumSCS, int ioProcNumAll, int scsMyId, MPI_Comm scsComm);

protected:
    //
    // Pointers to amr,amrlevel.
    //
    amrex::Amr*             parent;
    amrex::Array<amrex::AmrLevel*> LevelData;
    //
    // Pointers to grad_phi at previous and current time
    //
    amrex::Array< amrex::Array<std::unique_ptr<amrex::MultiFab> > > grad_phi_curr;
    amrex::Array< amrex::Array<std::unique_ptr<amrex::MultiFab> > > grad_phi_prev;

    amrex::Array<std::unique_ptr<amrex::FluxRegister> > phi_flux_reg;
    //
    // amrex::BoxArray at each level
    //
    const amrex::Array<amrex::BoxArray>& grids;
    //
    const amrex::Array<amrex::DistributionMapping>& dmap;
    //
    // Resnorm at each level
    //
    amrex::Array<Real> level_solver_resnorm;

    int density;
    int finest_level;
    int finest_level_allocated;

    int mg_bc[2*BL_SPACEDIM];

    amrex::BCRec* phys_bc;

    static int verbose;
    static int show_timings;
    static int no_sync;
    static int no_composite;
    static int dirichlet_bcs;
    static int  monopole_bcs;
    static int  solve_with_cpp;
    static int solve_with_hpgmg;
    static Real mass_offset;
    static Real sl_tol;
    static Real ml_tol;
    static Real delta_tol;
    static std::string gravity_type;
    static int stencil_type;

    void fill_ec_grow(int level, const amrex::Array<amrex::MultiFab*>& ecF,
                      const amrex::Array<amrex::MultiFab*>& ecC) const;

    void AddParticlesToRhs(int level, amrex::MultiFab& Rhs, int ngrow);
    void AddGhostParticlesToRhs(int level, amrex::MultiFab& Rhs);
    void AddVirtualParticlesToRhs(int level, amrex::MultiFab& Rhs, int ngrow);

    void AddParticlesToRhs(int base_level, int finest_level, const amrex::Array<amrex::MultiFab*>& Rhs_particles);
    void AddGhostParticlesToRhs(int level, const amrex::Array<amrex::MultiFab*>& Rhs_particles);
    void AddVirtualParticlesToRhs(int finest_level, const amrex::Array<amrex::MultiFab*>& Rhs_particles);

    void CorrectRhsUsingOffset(int level, amrex::MultiFab& Rhs);
};
#endif

