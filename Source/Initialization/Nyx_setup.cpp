
#include <AMReX_LevelBld.H>
#include <AMReX_buildInfo.H>

#include <Nyx.H>
#ifndef NO_HYDRO
#include <Derive.H>
#endif
#include <ParticleDerive.H>
#include <bc_fill.H>

#ifdef HEATCOOL
#include <atomic_rates.H>
#endif

#include <Forcing.H>

using namespace amrex;
using std::string;

static Box the_same_box(const Box& b)
{
    return b;
}

static Box grow_box_by_one(const Box& b)
{
    return amrex::grow(b, 1);
}

typedef StateDescriptor::BndryFunc BndryFunc;

//
// Components are:
//  Interior, Inflow, Outflow,  Symmetry,     SlipWall,     NoSlipWall
//
static int scalar_bc[] =
{
    INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN
};

static int norm_vel_bc[] =
{
    INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_ODD, REFLECT_ODD, REFLECT_ODD
};

static int tang_vel_bc[] =
{
    INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN
};

static
void
set_scalar_bc(BCRec& bc, const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    for (int i = 0; i < AMREX_SPACEDIM; i++)
    {
        bc.setLo(i, scalar_bc[lo_bc[i]]);
        bc.setHi(i, scalar_bc[hi_bc[i]]);
    }
}

static
void
set_x_vel_bc(BCRec& bc, const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    bc.setLo(0, norm_vel_bc[lo_bc[0]]);
    bc.setHi(0, norm_vel_bc[hi_bc[0]]);
    bc.setLo(1, tang_vel_bc[lo_bc[1]]);
    bc.setHi(1, tang_vel_bc[hi_bc[1]]);
    bc.setLo(2, tang_vel_bc[lo_bc[2]]);
    bc.setHi(2, tang_vel_bc[hi_bc[2]]);
}

static
void
set_y_vel_bc(BCRec& bc, const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    bc.setLo(0, tang_vel_bc[lo_bc[0]]);
    bc.setHi(0, tang_vel_bc[hi_bc[0]]);
    bc.setLo(1, norm_vel_bc[lo_bc[1]]);
    bc.setHi(1, norm_vel_bc[hi_bc[1]]);
    bc.setLo(2, tang_vel_bc[lo_bc[2]]);
    bc.setHi(2, tang_vel_bc[hi_bc[2]]);
}

static
void
set_z_vel_bc(BCRec& bc, const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    bc.setLo(0, tang_vel_bc[lo_bc[0]]);
    bc.setHi(0, tang_vel_bc[hi_bc[0]]);
    bc.setLo(1, tang_vel_bc[lo_bc[1]]);
    bc.setHi(1, tang_vel_bc[hi_bc[1]]);
    bc.setLo(2, norm_vel_bc[lo_bc[2]]);
    bc.setHi(2, norm_vel_bc[hi_bc[2]]);
}

void
Nyx::variable_setup()
{

  // initialize the start time for our CPU-time tracker
  startCPUTime = ParallelDescriptor::second();

    BL_ASSERT(desc_lst.size() == 0);

    if (ParallelDescriptor::IOProcessor()) {
          const char* amrex_hash  = amrex::buildInfoGetGitHash(2);
          std::cout << "\n" << "AMReX git describe: " << amrex_hash << "\n";
          const char* nyx_hash  = amrex::buildInfoGetGitHash(1);
          std::cout << "\n" << "Nyx git describe:   " << nyx_hash << "\n";
    }

    // Get options, set phys_bc
    read_params();

#ifdef NO_HYDRO
    if (do_dm_particles)
        no_hydro_setup();
#ifdef HEATCOOL
    else
        heatcool_setup();
#endif
#else
    if (do_hydro == 1)
    {
       hydro_setup();
    }
    else if (do_dm_particles)
    {
       no_hydro_setup();
    }
#ifdef HEATCOOL
    else
        heatcool_setup();
#endif
#endif

    //
    // DEFINE ERROR ESTIMATION QUANTITIES
    //
    error_setup();
}

#ifndef NO_HYDRO
#ifdef HEATCOOL
void
Nyx::heatcool_setup ()
{
    ParmParse pp_nyx("nyx");
    std::string file_in;
    pp_nyx.query("uvb_rates_file", file_in);
    amrex::Real mean_rhob = comoving_OmB * 3.e0*(comoving_h*100.e0)*(comoving_h*100.e0) / (8.e0*M_PI*Gconst);
    tabulate_rates(file_in, mean_rhob);
    amrex::Gpu::streamSynchronize();
}
#endif

void
Nyx::hydro_setup()
{
    //
    // Set number of state variables and pointers to components
    //
    int cnt = 6;

    // Note that we must set NDIAG_C before we call set_method_params 
    int NDIAG_C;
    if (inhomo_reion > 0)
    {
        NDIAG_C  = 3;
        Zhi_comp = 2;
    } else {
        NDIAG_C  = 2;
    }

#ifdef CONST_SPECIES
    NumSpec = 0;
#else
    // Get the number of species from the network model.
    NumSpec = 2;
    if (NumSpec > 0)
        cnt += NumSpec;
#endif

    NUM_STATE = cnt;


#ifdef HEATCOOL
    heatcool_setup();
#endif

    Interpolater* interp = &cell_cons_interp;

    // Note that the default is state_data_extrap = false,
    // store_in_checkpoint = true.  We only need to put these in
    // explicitly if we want to do something different,
    // like not store the state data in a checkpoint directory
    bool state_data_extrap = false;
    bool store_in_checkpoint;

    store_in_checkpoint = true;
    desc_lst.addDescriptor(State_Type, IndexType::TheCellType(),
                           StateDescriptor::Point, nghost_state, NUM_STATE, interp,
                           state_data_extrap, store_in_checkpoint);

    // This has two components: Temperature and Ne
    desc_lst.addDescriptor(DiagEOS_Type, IndexType::TheCellType(),
                           StateDescriptor::Point, nghost_state, NDIAG_C, interp,
                           state_data_extrap, store_in_checkpoint);

#ifdef SDC
    // This only has one component -- the update to rho_e from reactions
    store_in_checkpoint = true;
    desc_lst.addDescriptor(SDC_IR_Type, IndexType::TheCellType(),
                           StateDescriptor::Point, 1, 1, interp,
                           state_data_extrap, store_in_checkpoint);
#endif

    store_in_checkpoint = true;
    desc_lst.addDescriptor(PhiGrav_Type, IndexType::TheCellType(),
                           StateDescriptor::Point, 1, 1,
                           &cell_cons_interp, state_data_extrap,
                           store_in_checkpoint);

    store_in_checkpoint = false;
    desc_lst.addDescriptor(Gravity_Type, IndexType::TheCellType(),
                           StateDescriptor::Point, 1, AMREX_SPACEDIM,
                           &cell_cons_interp, state_data_extrap,
                           store_in_checkpoint);

    Vector<BCRec> bcs(NUM_STATE);
    Vector<std::string> name(NUM_STATE);

    BCRec bc;
    cnt = 0;
    set_scalar_bc(bc, phys_bc);  bcs[cnt] = bc;  name[cnt] = "density";
    cnt++;
    set_x_vel_bc(bc, phys_bc);  bcs[cnt] = bc;  name[cnt] = "xmom";
    cnt++;
    set_y_vel_bc(bc, phys_bc);  bcs[cnt] = bc;  name[cnt] = "ymom";
    cnt++;
    set_z_vel_bc(bc, phys_bc);  bcs[cnt] = bc;  name[cnt] = "zmom";
    cnt++;
    set_scalar_bc(bc, phys_bc);  bcs[cnt] = bc;  name[cnt] = "rho_E";
    cnt++;
    set_scalar_bc(bc, phys_bc);  bcs[cnt] = bc;  name[cnt] = "rho_e";

    // Get the species names from the network model.
    Vector<std::string> spec_names(2);
    spec_names[0] = "H";
    spec_names[1] = "He";

    if (ParallelDescriptor::IOProcessor())
    {
        std::cout << NumSpec << " Species: ";
        for (int i = 0; i < NumSpec; i++)
           std::cout << spec_names[i] << ' ' << ' ';
        std::cout << '\n';
    }


#ifndef CONST_SPECIES
    for (int i = 0; i < NumSpec; ++i)
    {
        cnt++;
        set_scalar_bc(bc,phys_bc);
        bcs[cnt] = bc;
        name[cnt] = "rho_" + spec_names[i];
     }
#endif

    StateDescriptor::BndryFunc bndryfunc(nyx_bcfill);
    bndryfunc.setRunOnGPU(true);  // I promise the bc function will launch gpu kernels.

    desc_lst.setComponent(State_Type,
                          Density_comp,
                          name,
                          bcs,
                          bndryfunc);

    set_scalar_bc(bc, phys_bc);
    desc_lst.setComponent(DiagEOS_Type, 0, "Temp", bc,
                          bndryfunc);
    desc_lst.setComponent(DiagEOS_Type, 1, "Ne", bc,
                          bndryfunc);

    if (inhomo_reion > 0) {
       desc_lst.setComponent(DiagEOS_Type, 2, "Z_HI", bc,
                             bndryfunc);
    }

#ifdef SDC
    set_scalar_bc(bc, phys_bc);
    desc_lst.setComponent(SDC_IR_Type, 0, "I_R", bc,
                          bndryfunc);
#endif

    if (do_grav)
    {
        set_scalar_bc(bc, phys_bc);
        desc_lst.setComponent(PhiGrav_Type, 0, "phi_grav", bc,
                              bndryfunc);

        set_x_vel_bc(bc, phys_bc);
        desc_lst.setComponent(Gravity_Type, 0, "grav_x", bc,
                             bndryfunc);

        set_y_vel_bc(bc, phys_bc);
        desc_lst.setComponent(Gravity_Type, 1, "grav_y", bc,
                              bndryfunc);

        set_z_vel_bc(bc, phys_bc);
        desc_lst.setComponent(Gravity_Type, 2, "grav_z", bc,
                              bndryfunc);
    }

    //
    // DEFINE DERIVED QUANTITIES
    //

    //
    // Density * Volume
    //
    derive_lst.add("denvol", IndexType::TheCellType(), 1,
                   derdenvol, grow_box_by_one);
    derive_lst.addComponent("denvol", desc_lst, State_Type, Density_comp, 1);

    //
    // Density / (avg_gas_density * 8^(level+1))
    //
    derive_lst.add("overden", IndexType::TheCellType(), 1,
                   deroverden, grow_box_by_one);
    derive_lst.addComponent("overden", desc_lst, State_Type, Density_comp, 1);

    //
    // Cell in center refined area
    //
    derive_lst.add("overdenzoom", IndexType::TheCellType(), 1,
                   deroverden, grow_box_by_one);
    derive_lst.addComponent("overdenzoom", desc_lst, State_Type, Density_comp, 1);

    //
    // Pressure
    //
    derive_lst.add("pressure", IndexType::TheCellType(), 1,
                   derpres, the_same_box);
    derive_lst.addComponent("pressure", desc_lst, State_Type, Density_comp,
                            NUM_STATE);

    //
    // Kinetic energy
    //
    derive_lst.add("kineng", IndexType::TheCellType(), 1,
                   derkineng, the_same_box);
    derive_lst.addComponent("kineng", desc_lst, State_Type, Density_comp, 1);
    derive_lst.addComponent("kineng", desc_lst, State_Type, Xmom_comp, AMREX_SPACEDIM);

    //
    // Sound speed (c)
    //
    derive_lst.add("soundspeed", IndexType::TheCellType(), 1,
                   dersoundspeed, the_same_box);
    derive_lst.addComponent("soundspeed", desc_lst, State_Type, Density_comp,
                            NUM_STATE);

    //
    // Mach number(M)
    //
    derive_lst.add("MachNumber", IndexType::TheCellType(), 1,
                   dermachnumber, the_same_box);
    derive_lst.addComponent("MachNumber", desc_lst, State_Type, Density_comp,
                            NUM_STATE);

    //
    // Gravitational forcing
    //
    //if (do_grav)  
    //{
    //derive_lst.add("rhog",IndexType::TheCellType(),1,
    //               rhog,the_same_box);
    //derive_lst.addComponent("rhog",desc_lst,State_Type,Density_comp,1);
    //derive_lst.addComponent("rhog",desc_lst,Gravity_Type,0,AMREX_SPACEDIM);
    //}

    //
    // Div(u)
    //
    derive_lst.add("divu", IndexType::TheCellType(), 1,
                   derdivu, grow_box_by_one);
    derive_lst.addComponent("divu", desc_lst, State_Type, Density_comp, 1);
    derive_lst.addComponent("divu", desc_lst, State_Type, Xmom_comp, AMREX_SPACEDIM);

    //
    // Internal energy as derived from rho*E, part of the state
    //
    derive_lst.add("eint_E", IndexType::TheCellType(), 1,
                   dereint1, the_same_box);
    derive_lst.addComponent("eint_E", desc_lst, State_Type, Density_comp, NUM_STATE);

    //
    // Internal energy as derived from rho*e, part of the state
    //
    derive_lst.add("eint_e", IndexType::TheCellType(), 1,
                   dereint2, the_same_box);
    derive_lst.addComponent("eint_e", desc_lst, State_Type, Density_comp, NUM_STATE);

    //
    // Log(density)
    //
    derive_lst.add("logden", IndexType::TheCellType(), 1,
                   derlogden, the_same_box);
    derive_lst.addComponent("logden", desc_lst, State_Type, Density_comp, 1);

    derive_lst.add("StateErr", IndexType::TheCellType(), 3,
                   derstate,
                   grow_box_by_one);
    derive_lst.addComponent("StateErr", desc_lst,   State_Type, Density_comp, 1);
    derive_lst.addComponent("StateErr", desc_lst, DiagEOS_Type, Temp_comp, 1);
#ifndef CONST_SPECIES
    derive_lst.addComponent("StateErr", desc_lst, State_Type, FirstSpec_comp, 1);
#endif

    //
    // X from rhoX
    //
#ifndef CONST_SPECIES
    for (int i = 0; i < NumSpec; i++)
    {
        string spec_string = "X(" + spec_names[i] + ")";

        derive_lst.add(spec_string, IndexType::TheCellType(), 1,
                       derspec, the_same_box);
        derive_lst.addComponent(spec_string, desc_lst, State_Type, Density_comp, 1);
        derive_lst.addComponent(spec_string, desc_lst, State_Type,
                                FirstSpec_comp + i, 1);
    }
#endif

    derive_lst.add("forcex", IndexType::TheCellType(), 1,
                   derforcex, the_same_box);
    derive_lst.addComponent("forcex", desc_lst, State_Type, Density_comp, 1);

    derive_lst.add("forcey", IndexType::TheCellType(), 1,
                   derforcey, the_same_box);
    derive_lst.addComponent("forcey", desc_lst, State_Type, Density_comp, 1);

    derive_lst.add("forcez", IndexType::TheCellType(), 1,
                   derforcez, the_same_box);
    derive_lst.addComponent("forcez", desc_lst, State_Type, Density_comp, 1);

    //
    // Velocities
    //
    derive_lst.add("x_velocity", IndexType::TheCellType(), 1,
                   dervel, the_same_box);
    derive_lst.addComponent("x_velocity", desc_lst, State_Type, Density_comp, 1);
    derive_lst.addComponent("x_velocity", desc_lst, State_Type, Xmom_comp, 1);

    derive_lst.add("y_velocity", IndexType::TheCellType(), 1,
                   dervel, the_same_box);
    derive_lst.addComponent("y_velocity", desc_lst, State_Type, Density_comp, 1);
    derive_lst.addComponent("y_velocity", desc_lst, State_Type, Ymom_comp, 1);

    derive_lst.add("z_velocity", IndexType::TheCellType(), 1,
                   dervel, the_same_box);
    derive_lst.addComponent("z_velocity", desc_lst, State_Type, Density_comp, 1);
    derive_lst.addComponent("z_velocity", desc_lst, State_Type, Zmom_comp, 1);

    //
    // Magnitude of velocity.
    //
    derive_lst.add("magvel", IndexType::TheCellType(), 1,
                   dermagvel, the_same_box);
    derive_lst.addComponent("magvel", desc_lst, State_Type, Density_comp, 1);
    derive_lst.addComponent("magvel", desc_lst, State_Type, Xmom_comp, AMREX_SPACEDIM);

    //
    // Magnitude of vorticity.
    //
    derive_lst.add("magvort",IndexType::TheCellType(),1,
                   dermagvort,grow_box_by_one);
    // Here we exploit the fact that Xmom = Density + 1
    //   in order to use the correct interpolation.
    if (Xmom_comp != Density_comp+1)
       amrex::Error("We are assuming Xmom = Density + 1 in Nyx_setup.cpp");
    derive_lst.addComponent("magvort",desc_lst,State_Type,Density_comp,AMREX_SPACEDIM+1);

    //
    // Magnitude of momentum.
    //
    derive_lst.add("magmom", IndexType::TheCellType(), 1,
                   dermagmom, the_same_box);
    derive_lst.addComponent("magmom", desc_lst, State_Type, Xmom_comp, AMREX_SPACEDIM);

    if (do_grav)  
    {
         derive_lst.add("maggrav", IndexType::TheCellType(), 1,
                        dermaggrav,
                        the_same_box);
         derive_lst.addComponent("maggrav", desc_lst, Gravity_Type, 0, AMREX_SPACEDIM);
    }

    //
    // We want a derived type that corresponds to the number of particles in
    // each cell. We only intend to use it in plotfiles for debugging purposes.
    // We'll just use the DERNULL since don't do anything in fortran for now.
    // We'll actually set the values in `writePlotFile()`.
    //
    derive_lst.add("particle_count", IndexType::TheCellType(), 1,
                   dernull, the_same_box);
    derive_lst.addComponent("particle_count", desc_lst, State_Type, Density_comp, 1);

    derive_lst.add("particle_mass_density", IndexType::TheCellType(), 1,
                   dernull, grow_box_by_one);
    derive_lst.addComponent("particle_mass_density", desc_lst, State_Type,
                            Density_comp, 1);

    derive_lst.add("particle_x_velocity", IndexType::TheCellType(), 1,
                   dernull, grow_box_by_one);
    derive_lst.addComponent("particle_x_velocity", desc_lst, State_Type,
                            Density_comp, 1);
    derive_lst.add("particle_y_velocity", IndexType::TheCellType(), 1,
                   dernull, grow_box_by_one);
    derive_lst.addComponent("particle_y_velocity", desc_lst, State_Type,
                            Density_comp, 1);
    derive_lst.add("particle_z_velocity", IndexType::TheCellType(), 1,
                   dernull, grow_box_by_one);
    derive_lst.addComponent("particle_z_velocity", desc_lst, State_Type,
                            Density_comp, 1);

#ifdef AGN
    derive_lst.add("agn_particle_count", IndexType::TheCellType(), 1,
                   dernull, the_same_box);
    derive_lst.addComponent("agn_particle_count", desc_lst, State_Type, Density_comp, 1);

    derive_lst.add("agn_mass_density", IndexType::TheCellType(), 1,
                   dernull, grow_box_by_one);
    derive_lst.addComponent("agn_mass_density", desc_lst, State_Type,
                            Density_comp, 1);
#endif

#ifdef NEUTRINO_PARTICLES
    derive_lst.add("neutrino_particle_count", IndexType::TheCellType(), 1,
                   dernull, the_same_box);
    derive_lst.addComponent("neutrino_particle_count", desc_lst, State_Type, Density_comp, 1);

    derive_lst.add("neutrino_mass_density", IndexType::TheCellType(), 1,
                   dernull, grow_box_by_one);
    derive_lst.addComponent("neutrino_mass_density", desc_lst, State_Type,
                            Density_comp, 1);

    derive_lst.add("neutrino_x_velocity", IndexType::TheCellType(), 1,
                   dernull, grow_box_by_one);
    derive_lst.addComponent("neutrino_x_velocity", desc_lst, State_Type,
                            Density_comp, 1);
    derive_lst.add("neutrino_y_velocity", IndexType::TheCellType(), 1,
                   dernull, grow_box_by_one);
    derive_lst.addComponent("neutrino_y_velocity", desc_lst, State_Type,
                            Density_comp, 1);
    derive_lst.add("neutrino_z_velocity", IndexType::TheCellType(), 1,
                   dernull, grow_box_by_one);
    derive_lst.addComponent("neutrino_z_velocity", desc_lst, State_Type,
                            Density_comp, 1);
#endif

    derive_lst.add("total_particle_count", IndexType::TheCellType(), 1,
                   dernull, the_same_box);
    derive_lst.addComponent("total_particle_count", desc_lst, State_Type,
                            Density_comp, 1);

    derive_lst.add("total_density", IndexType::TheCellType(), 1,
                   dernull, grow_box_by_one);
    derive_lst.addComponent("total_density", desc_lst, State_Type,
                            Density_comp, 1);
    derive_lst.add("Rank", IndexType::TheCellType(), 1,
                   dernull, grow_box_by_one);

#ifndef CONST_SPECIES
    for (int i = 0; i < NumSpec; i++)
    {
        derive_lst.add(spec_names[i], IndexType::TheCellType(), 1,
                       derspec, the_same_box);
        derive_lst.addComponent(spec_names[i], desc_lst, State_Type, Density_comp, 1);
        derive_lst.addComponent(spec_names[i], desc_lst, State_Type,
                                FirstSpec_comp + i, 1);
    }
#endif
}
#endif

void
Nyx::no_hydro_setup()
{
    NUM_STATE = 1;

#ifdef HEATCOOL
    heatcool_setup();
#endif

    // Note that the default is state_data_extrap = false,
    // store_in_checkpoint = true.  We only need to put these in
    // explicitly if we want to do something different,
    // like not store the state data in a checkpoint directory
    bool state_data_extrap = false;
    bool store_in_checkpoint;

    BCRec bc;

    StateDescriptor::BndryFunc bndryfunc(nyx_bcfill);
    bndryfunc.setRunOnGPU(true);  // I promise the bc function will launch gpu kernels.

#ifndef NO_HYDRO
    // We have to create these anyway because the StateTypes are defined at compile time.
    // However, we define them with only one component each because we don't actually use them.
    store_in_checkpoint = false;
    // This has only one dummy components
    desc_lst.addDescriptor(State_Type, IndexType::TheCellType(),
                           StateDescriptor::Point, 1, 1, &cell_cons_interp,
                           state_data_extrap, store_in_checkpoint);

    set_scalar_bc(bc, phys_bc);
    desc_lst.setComponent(State_Type, 0, "density", bc,
                          bndryfunc);

    // This has only one dummy components
    store_in_checkpoint = false;
    desc_lst.addDescriptor(DiagEOS_Type, IndexType::TheCellType(),
                           StateDescriptor::Point, 1, 1, &cell_cons_interp,
                           state_data_extrap, store_in_checkpoint);

    set_scalar_bc(bc, phys_bc);
    desc_lst.setComponent(DiagEOS_Type, 0, "Temp", bc,
                          bndryfunc);
#endif

    store_in_checkpoint = true;
    desc_lst.addDescriptor(PhiGrav_Type, IndexType::TheCellType(),
                           StateDescriptor::Point, 1, 1,
                           &cell_cons_interp, state_data_extrap,
                           store_in_checkpoint);

    store_in_checkpoint = false;
    desc_lst.addDescriptor(Gravity_Type, IndexType::TheCellType(),
                           StateDescriptor::Point, 1, AMREX_SPACEDIM,
                           &cell_cons_interp, state_data_extrap,
                           store_in_checkpoint);

    if (do_grav)
    {
       set_scalar_bc(bc, phys_bc);
       desc_lst.setComponent(PhiGrav_Type, 0, "phi_grav", bc,
                             bndryfunc);
       set_x_vel_bc(bc, phys_bc);
       desc_lst.setComponent(Gravity_Type, 0, "grav_x", bc,
                             bndryfunc);
       set_y_vel_bc(bc, phys_bc);
       desc_lst.setComponent(Gravity_Type, 1, "grav_y", bc,
                             bndryfunc);
       set_z_vel_bc(bc, phys_bc);
       desc_lst.setComponent(Gravity_Type, 2, "grav_z", bc,
                             bndryfunc);

       derive_lst.add("maggrav", IndexType::TheCellType(), 1,
                      dermaggrav,
                      the_same_box);
       derive_lst.addComponent("maggrav", desc_lst, Gravity_Type, 0, AMREX_SPACEDIM);
    }

    // We want a derived type that corresponds to the number of particles in
    // each cell. We only intend to use it in plotfiles for debugging purposes.
    // We'll just use the DERNULL since don't do anything in fortran for now.
    // We'll actually set the values in `writePlotFile()`.
    //
    derive_lst.add("particle_count", IndexType::TheCellType(), 1,
                   dernull, the_same_box);
    derive_lst.addComponent("particle_count", desc_lst, PhiGrav_Type, 0, 1);

    derive_lst.add("particle_mass_density", IndexType::TheCellType(), 1,
                   dernull, grow_box_by_one);
    derive_lst.addComponent("particle_mass_density", desc_lst, PhiGrav_Type, 0, 1);
#ifndef NO_HYDRO
    derive_lst.add("particle_x_velocity", IndexType::TheCellType(), 1,
                   dernull, grow_box_by_one);
    derive_lst.addComponent("particle_x_velocity", desc_lst, State_Type,
                            Density_comp, 1);
    derive_lst.add("particle_y_velocity", IndexType::TheCellType(), 1,
                   dernull, grow_box_by_one);
    derive_lst.addComponent("particle_y_velocity", desc_lst, State_Type,
                            Density_comp, 1);
    derive_lst.add("particle_z_velocity", IndexType::TheCellType(), 1,
                   dernull, grow_box_by_one);
    derive_lst.addComponent("particle_z_velocity", desc_lst, State_Type,
                            Density_comp, 1);
    //
    // Density * Volume
    //
    derive_lst.add("denvol", IndexType::TheCellType(), 1,
                   derdenvol, grow_box_by_one);
    derive_lst.addComponent("denvol", desc_lst, State_Type, Density_comp, 1);

    //
    // Density / (avg_gas_density * 8^(level+1))
    //
    derive_lst.add("overden", IndexType::TheCellType(), 1,
                   deroverden, grow_box_by_one);
    derive_lst.addComponent("overden", desc_lst, State_Type, Density_comp, 1);

    //
    // Cell in center refined area
    //
    derive_lst.add("overdenzoom", IndexType::TheCellType(), 1,
                   deroverden, grow_box_by_one);
    derive_lst.addComponent("overdenzoom", desc_lst, State_Type, Density_comp, 1);

#else
    derive_lst.add("particle_x_velocity", IndexType::TheCellType(), 1,
                   dernull, grow_box_by_one);
    derive_lst.addComponent("particle_x_velocity", desc_lst, PhiGrav_Type,
                            0, 1);
    derive_lst.add("particle_y_velocity", IndexType::TheCellType(), 1,
                   dernull, grow_box_by_one);
    derive_lst.addComponent("particle_y_velocity", desc_lst, PhiGrav_Type,
                            0, 1);
    derive_lst.add("particle_z_velocity", IndexType::TheCellType(), 1,
                   dernull, grow_box_by_one);
    derive_lst.addComponent("particle_z_velocity", desc_lst, PhiGrav_Type,
                            0, 1);
    //
    // Density * Volume
    //
    derive_lst.add("denvol", IndexType::TheCellType(), 1,
                   derdenvol, grow_box_by_one);
    derive_lst.addComponent("denvol", desc_lst, PhiGrav_Type, 0, 1);

    //
    // Density / (avg_gas_density * 8^(level+1))
    //
    derive_lst.add("overden", IndexType::TheCellType(), 1,
                   deroverden, grow_box_by_one);
    derive_lst.addComponent("overden", desc_lst, PhiGrav_Type, 0, 1);

    //
    // Cell in center refined area
    //
    derive_lst.add("overdenzoom", IndexType::TheCellType(), 1,
                   deroverden, grow_box_by_one);
    derive_lst.addComponent("overdenzoom", desc_lst, PhiGrav_Type, 0, 1);

#endif
    derive_lst.add("total_particle_count", IndexType::TheCellType(), 1,
                   dernull, grow_box_by_one);
    derive_lst.addComponent("total_particle_count", desc_lst, PhiGrav_Type, 0, 1);

    derive_lst.add("total_density", IndexType::TheCellType(), 1,
                   dernull, grow_box_by_one);
    derive_lst.addComponent("total_density", desc_lst, PhiGrav_Type, 0, 1);

#ifdef AGN
    derive_lst.add("agn_mass_density", IndexType::TheCellType(), 1,
                   dernull, grow_box_by_one);
    derive_lst.addComponent("agn_mass_density", desc_lst, Gravity_Type, 0, 1);
#endif
#ifdef NEUTRINO_PARTICLES
#ifndef NO_HYDRO
    derive_lst.add("neutrino_particle_count", IndexType::TheCellType(), 1,
                   dernull, grow_box_by_one);
    derive_lst.addComponent("neutrino_particle_count", desc_lst, State_Type, Density_comp, 1);

    derive_lst.add("neutrino_mass_density", IndexType::TheCellType(), 1,
                   dernull, grow_box_by_one);
    derive_lst.addComponent("neutrino_mass_density", desc_lst, State_Type,
                            Density_comp, 1);
    derive_lst.add("neutrino_x_velocity", IndexType::TheCellType(), 1,
                   dernull, grow_box_by_one);
    derive_lst.addComponent("neutrino_x_velocity", desc_lst, State_Type,
                            Density_comp, 1);
    derive_lst.add("neutrino_y_velocity", IndexType::TheCellType(), 1,
                   dernull, grow_box_by_one);
    derive_lst.addComponent("neutrino_y_velocity", desc_lst, State_Type,
                            Density_comp, 1);
    derive_lst.add("neutrino_z_velocity", IndexType::TheCellType(), 1,
                   dernull, grow_box_by_one);
    derive_lst.addComponent("neutrino_z_velocity", desc_lst, State_Type,
                            Density_comp, 1);
#else
    derive_lst.add("neutrino_particle_count", IndexType::TheCellType(), 1,
                   dernull, grow_box_by_one);
    derive_lst.addComponent("neutrino_particle_count", desc_lst, PhiGrav_Type, 0, 1);

    derive_lst.add("neutrino_mass_density", IndexType::TheCellType(), 1,
                   dernull, grow_box_by_one);
    derive_lst.addComponent("neutrino_mass_density", desc_lst, PhiGrav_Type,
                            0, 1);
    derive_lst.add("neutrino_x_velocity", IndexType::TheCellType(), 1,
                   dernull, grow_box_by_one);
    derive_lst.addComponent("neutrino_x_velocity", desc_lst, PhiGrav_Type,
                            0, 1);
    derive_lst.add("neutrino_y_velocity", IndexType::TheCellType(), 1,
                   dernull, grow_box_by_one);
    derive_lst.addComponent("neutrino_y_velocity", desc_lst, PhiGrav_Type,
                            0, 1);
    derive_lst.add("neutrino_z_velocity", IndexType::TheCellType(), 1,
                   dernull, grow_box_by_one);
    derive_lst.addComponent("neutrino_z_velocity", desc_lst, PhiGrav_Type,
                            0, 1);
#endif
#endif
}
