#include <winstd.H>

#include "LevelBld.H"

#include "Nyx.H"
#include "Nyx_F.H"
#include "Derive_F.H"

using std::string;

static Box the_same_box(const Box& b)
{
    return b;
}

static Box grow_box_by_one(const Box& b)
{
    return BoxLib::grow(b, 1);
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
    for (int i = 0; i < BL_SPACEDIM; i++)
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
    BL_ASSERT(desc_lst.size() == 0);

    // Initialize the network
    network_init();

    // Get options, set phys_bc
    read_params();

    if (do_hydro == 1) 
    {
       hydro_setup();
    }
    else
    {
       no_hydro_setup();
    }

    //
    // DEFINE ERROR ESTIMATION QUANTITIES
    //
    error_setup();
}

void
Nyx::hydro_setup()
{
    //
    // Set number of state variables and pointers to components
    //
    int use_sgs = 0;
    int use_ssfm = 0;
    Density = 0;
    Xmom    = 1;
    Ymom    = 2;
    Zmom    = 3;
    Eden    = 4;
    Eint    = 5;
    int cnt = 6;
#ifdef SGS
    Esgs = cnt++;
    use_sgs = 1;
#endif
    Temp = cnt++;

    NumAdv = 1;
    if (NumAdv > 0)
    {
        FirstAdv = cnt;
        cnt += NumAdv;
    }

#ifdef SSFM
    if(NumAdv == 0){
      FirstAdv = cnt;
    }
    NumAdv += 2;
    use_ssfm = 1;
    cnt += 2;
#endif


    int dm = BL_SPACEDIM;

    // Get the number of species from the network model.
    BL_FORT_PROC_CALL(GET_NUM_SPEC, get_num_spec)(&NumSpec);

    if (NumSpec > 0)
    {
        FirstSpec = cnt;
        cnt += NumSpec;
    }

    // Get the number of species from the network model.
    BL_FORT_PROC_CALL(GET_NUM_AUX, get_num_aux)(&NumAux);

    int FirstAux = -1;
    if (NumAux > 0)
    {
        FirstAux = cnt;
        cnt += NumAux;
    }

    NUM_STATE = cnt;

    // Define NUM_GROW from the f90 module.
    BL_FORT_PROC_CALL(GET_METHOD_PARAMS, get_method_params)(&NUM_GROW);

    BL_FORT_PROC_CALL(SET_METHOD_PARAMS, set_method_params)
        (dm, NumAdv, do_hydro, ppm_type, use_colglaz, 
         gamma, normalize_species ,use_sgs, use_ssfm);

    int coord_type = Geometry::Coord();
    BL_FORT_PROC_CALL(SET_PROBLEM_PARAMS, set_problem_params)
         (dm, phys_bc.lo(), phys_bc.hi(), Outflow, Symmetry, coord_type);

    Interpolater* interp = &cell_cons_interp;

    // Note that the default is state_data_extrap = false,
    // store_in_checkpoint = true.  We only need to put these in
    // explicitly if we want to do something different,
    // like not store the state data in a checkpoint directory
    bool state_data_extrap = false;
    bool store_in_checkpoint;

    store_in_checkpoint = true;
    desc_lst.addDescriptor(State_Type, IndexType::TheCellType(),
                           StateDescriptor::Point, 1, NUM_STATE, interp,
                           state_data_extrap, store_in_checkpoint);

#ifdef GRAVITY
    store_in_checkpoint = true;
    desc_lst.addDescriptor(PhiGrav_Type, IndexType::TheCellType(),
                           StateDescriptor::Point, 1, 1,
                           &cell_cons_interp, state_data_extrap,
                           store_in_checkpoint);

    store_in_checkpoint = false;
    desc_lst.addDescriptor(Gravity_Type, IndexType::TheCellType(),
                           StateDescriptor::Point, 1, BL_SPACEDIM,
                           &cell_cons_interp, state_data_extrap,
                           store_in_checkpoint);
#endif

#ifdef SGS
    // Component 0: prod_sgs
    // Component 1: diss_sgs
    // Component 2: turbulent forcing
    store_in_checkpoint = true;
    desc_lst.addDescriptor(SGS_Type, IndexType::TheCellType(),
                           StateDescriptor::Point, 1, 3, &cell_cons_interp,
                           state_data_extrap, store_in_checkpoint);
#endif

    Array<BCRec> bcs(NUM_STATE);
    Array<std::string> name(NUM_STATE);

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
#ifdef SGS
    cnt++;
    set_scalar_bc(bc, phys_bc);  bcs[cnt] = bc;  name[cnt] = "rho_K";
#endif
    cnt++;
    set_scalar_bc(bc, phys_bc);  bcs[cnt] = bc;  name[cnt] = "Temp";

#ifdef SSFM
    NumAdv -= 2;
#endif
    for (int i = 0; i < NumAdv; ++i)
    {
        cnt++;
        set_scalar_bc(bc, phys_bc);
        bcs[cnt]  = bc;
        name[cnt] = BoxLib::Concatenate("adv_", i, 1);
    }
#ifdef SSFM
    NumAdv += 2;
    cnt++;
    set_scalar_bc(bc, phys_bc); bcs[cnt] = bc; name[cnt] = "star_dens";
    cnt++;
    set_scalar_bc(bc, phys_bc); bcs[cnt] = bc; name[cnt] = "FB_res";
#endif
    // Get the species names from the network model.
    Array<std::string> spec_names(NumSpec);

    for (int i = 0; i < NumSpec; i++)
    {
        int len = 20;
        Array<int> int_spec_names(len);

        // This call return the actual length of each string in "len"
        BL_FORT_PROC_CALL(GET_SPEC_NAMES, get_spec_names)
            (int_spec_names.dataPtr(), &i, &len);

        for (int j = 0; j < len; j++)
            spec_names[i].push_back(int_spec_names[j]);
    }

    if (ParallelDescriptor::IOProcessor())
    {
        std::cout << NumSpec << " Species: ";
        for (int i = 0; i < NumSpec; i++)
           std::cout << spec_names[i] << ' ' << ' ';
        std::cout << '\n';
    }

    for (int i = 0; i < NumSpec; ++i)
    {
        cnt++;
        set_scalar_bc(bc,phys_bc);
        bcs[cnt] = bc;
        name[cnt] = "rho_" + spec_names[i];
    }

    // Get the auxiliary names from the network model.
    Array<std::string> aux_names(NumAux);

    for (int i = 0; i < NumAux; i++)
    {
        int len = 20;
        Array<int> int_aux_names(len);

        // This call return the actual length of each string in "len"
        BL_FORT_PROC_CALL(GET_AUX_NAMES, get_aux_names)
            (int_aux_names.dataPtr(), &i, &len);

        for (int j = 0; j < len; j++)
            aux_names[i].push_back(int_aux_names[j]);
    }

    if (ParallelDescriptor::IOProcessor())
    {
        std::cout << NumAux << " Auxiliary Variables: ";
        for (int i = 0; i < NumAux; i++)
           std::cout << aux_names[i] << ' ' << ' ';
        std::cout << '\n';
    }

    for (int i = 0; i < NumAux; ++i)
    {
        cnt++;
        set_scalar_bc(bc, phys_bc);
        bcs[cnt] = bc;
        name[cnt] = "rho_" + aux_names[i];
    }

    desc_lst.setComponent(State_Type, Density, name, bcs,
                          BndryFunc(BL_FORT_PROC_CALL(CA_DENFILL, ca_denfill),
                                    BL_FORT_PROC_CALL(CA_HYPFILL, ca_hypfill)));

#ifdef GRAVITY
    if (do_grav)
    {
        set_scalar_bc(bc, phys_bc);
        desc_lst.setComponent(PhiGrav_Type, 0, "phi_grav", bc,
                              BndryFunc(BL_FORT_PROC_CALL(CA_PHIFILL,
                                                          ca_phifill)));
        set_x_vel_bc(bc, phys_bc);
        desc_lst.setComponent(Gravity_Type, 0, "grav_x", bc,
                              BndryFunc(BL_FORT_PROC_CALL(CA_GRAVXFILL,
                                                          ca_gravxfill)));
       set_y_vel_bc(bc, phys_bc);
       desc_lst.setComponent(Gravity_Type, 1, "grav_y", bc,
                             BndryFunc(BL_FORT_PROC_CALL(CA_GRAVYFILL,
                                                         ca_gravyfill)));
       set_z_vel_bc(bc, phys_bc);
       desc_lst.setComponent(Gravity_Type, 2, "grav_z", bc,
                             BndryFunc(BL_FORT_PROC_CALL(CA_GRAVZFILL,
                                                         ca_gravzfill)));
    }
#endif

#ifdef SGS
    set_scalar_bc(bc, phys_bc);
    desc_lst.setComponent(SGS_Type, 0, "prod_sgs", bc,
                          BndryFunc(BL_FORT_PROC_CALL(CA_SGSFILL, ca_sgsfill)));
    desc_lst.setComponent(SGS_Type, 1, "diss_sgs", bc,
                          BndryFunc(BL_FORT_PROC_CALL(CA_SGSFILL, ca_sgsfill)));
    desc_lst.setComponent(SGS_Type, 2, "turb_src", bc,
                          BndryFunc(BL_FORT_PROC_CALL(CA_SGSFILL, ca_sgsfill)));
#endif

    //
    // DEFINE DERIVED QUANTITIES
    //
    // Pressure
    //
    derive_lst.add("pressure", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(CA_DERPRES, ca_derpres), the_same_box);
    derive_lst.addComponent("pressure", desc_lst, State_Type, Density,
                            NUM_STATE);

    //
    // Kinetic energy
    //
    derive_lst.add("kineng", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(CA_DERKINENG, ca_derkineng), the_same_box);
    derive_lst.addComponent("kineng", desc_lst, State_Type, Density, 1);
    derive_lst.addComponent("kineng", desc_lst, State_Type, Xmom, BL_SPACEDIM);

    //
    // Sound speed (c)
    //
    derive_lst.add("soundspeed", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(CA_DERSOUNDSPEED, ca_dersoundspeed),
                   the_same_box);
    derive_lst.addComponent("soundspeed", desc_lst, State_Type, Density,
                            NUM_STATE);

    //
    // Mach number(M)
    //
    derive_lst.add("MachNumber", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(CA_DERMACHNUMBER, ca_dermachnumber),
                   the_same_box);
    derive_lst.addComponent("MachNumber", desc_lst, State_Type, Density,
                            NUM_STATE);

    //
    // Gravitational forcing
    //
#ifdef GRAVITY
    //derive_lst.add("rhog",IndexType::TheCellType(),1,
    //               BL_FORT_PROC_CALL(CA_RHOG,ca_rhog),the_same_box);
    //derive_lst.addComponent("rhog",desc_lst,State_Type,Density,1);
    //derive_lst.addComponent("rhog",desc_lst,Gravity_Type,0,BL_SPACEDIM);
#endif

    //
    // Entropy (S)
    //
    derive_lst.add("entropy", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(CA_DERENTROPY, ca_derentropy),
                   the_same_box);
    derive_lst.addComponent("entropy", desc_lst, State_Type, Density,
                            NUM_STATE);

    //
    // Div(u)
    //
    derive_lst.add("divu", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(CA_DERDIVU, ca_derdivu), grow_box_by_one);
    derive_lst.addComponent("divu", desc_lst, State_Type, Density, 1);
    derive_lst.addComponent("divu", desc_lst, State_Type, Xmom, BL_SPACEDIM);

    //
    // Internal energy as derived from rho*E, part of the state
    //
    derive_lst.add("eint_E", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(CA_DEREINT1, ca_dereint1), the_same_box);
    derive_lst.addComponent("eint_E", desc_lst, State_Type, Density, NUM_STATE);

    //
    // Internal energy as derived from rho*e, part of the state
    //
    derive_lst.add("eint_e", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(CA_DEREINT2, ca_dereint2), the_same_box);
    derive_lst.addComponent("eint_e", desc_lst, State_Type, Density, NUM_STATE);

    //
    // Log(density)
    //
    derive_lst.add("logden", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(CA_DERLOGDEN, ca_derlogden), the_same_box);
    derive_lst.addComponent("logden", desc_lst, State_Type, Density, NUM_STATE);

    derive_lst.add("StateErr", IndexType::TheCellType(), 3,
                   BL_FORT_PROC_CALL(CA_DERSTATE, ca_derstate),
                   grow_box_by_one);
    derive_lst.addComponent("StateErr", desc_lst, State_Type, Density, 1);
    derive_lst.addComponent("StateErr", desc_lst, State_Type, Temp, 1);
    derive_lst.addComponent("StateErr", desc_lst, State_Type, FirstSpec, 1);

    //
    // X from rhoX
    //
    for (int i = 0; i < NumSpec; i++)
    {
        // @todo: smart string formatting
        string spec_string = "X(" + spec_names[i] + ")";

        derive_lst.add(spec_string, IndexType::TheCellType(), 1,
                       BL_FORT_PROC_CALL(CA_DERSPEC, ca_derspec), the_same_box);
        derive_lst.addComponent(spec_string, desc_lst, State_Type, Density, 1);
        derive_lst.addComponent(spec_string, desc_lst, State_Type,
                                FirstSpec + i, 1);
    }

    //
    // Forcing
    //
    derive_lst.add("forcex", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(CA_DERVEL, ca_dervel), the_same_box);
    derive_lst.addComponent("forcex", desc_lst, State_Type, Density, 1);
    derive_lst.addComponent("forcex", desc_lst, State_Type, Xmom, 1);

    derive_lst.add("forcey", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(CA_DERVEL, ca_dervel), the_same_box);
    derive_lst.addComponent("forcey", desc_lst, State_Type, Density, 1);
    derive_lst.addComponent("forcey", desc_lst, State_Type, Ymom, 1);

    derive_lst.add("forcez", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(CA_DERVEL, ca_dervel), the_same_box);
    derive_lst.addComponent("forcez", desc_lst, State_Type, Density, 1);
    derive_lst.addComponent("forcez", desc_lst, State_Type, Zmom, 1);

    //
    // Velocities
    //
    derive_lst.add("x_velocity", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(CA_DERVEL, ca_dervel), the_same_box);
    derive_lst.addComponent("x_velocity", desc_lst, State_Type, Density, 1);
    derive_lst.addComponent("x_velocity", desc_lst, State_Type, Xmom, 1);

    derive_lst.add("y_velocity", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(CA_DERVEL, ca_dervel), the_same_box);
    derive_lst.addComponent("y_velocity", desc_lst, State_Type, Density, 1);
    derive_lst.addComponent("y_velocity", desc_lst, State_Type, Ymom, 1);

    derive_lst.add("z_velocity", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(CA_DERVEL, ca_dervel), the_same_box);
    derive_lst.addComponent("z_velocity", desc_lst, State_Type, Density, 1);
    derive_lst.addComponent("z_velocity", desc_lst, State_Type, Zmom, 1);

#ifdef SGS
    derive_lst.add("K", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(CA_DERVEL, ca_dervel), the_same_box);
    derive_lst.addComponent("K", desc_lst, State_Type, Density, 1);
    derive_lst.addComponent("K", desc_lst, State_Type, Esgs, 1);
#endif

    //
    // Magnitude of velocity.
    //
    derive_lst.add("magvel", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(CA_DERMAGVEL, ca_dermagvel), the_same_box);
    derive_lst.addComponent("magvel", desc_lst, State_Type, Density, 1);
    derive_lst.addComponent("magvel", desc_lst, State_Type, Xmom, BL_SPACEDIM);

    //
    // Magnitude of vorticity.
    //
    derive_lst.add("magvort",IndexType::TheCellType(),1,
                   BL_FORT_PROC_CALL(CA_DERMAGVORT,ca_dermagvort),grow_box_by_one);
    // Here we exploit the fact that Xmom = Density + 1
    //   in order to use the correct interpolation.
    if (Xmom != Density+1)
       BoxLib::Error("We are assuming Xmom = Density + 1 in Nyx_setup.cpp");
    derive_lst.addComponent("magvort",desc_lst,State_Type,Density,BL_SPACEDIM+1);

    //
    // Magnitude of momentum.
    //
    derive_lst.add("magmom", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(CA_DERMAGMOM, ca_dermagmom), the_same_box);
    derive_lst.addComponent("magmom", desc_lst, State_Type, Xmom, BL_SPACEDIM);

#ifdef GRAVITY
    derive_lst.add("maggrav", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(CA_DERMAGGRAV, ca_dermaggrav),
                   the_same_box);
    derive_lst.addComponent("maggrav", desc_lst, Gravity_Type, 0, BL_SPACEDIM);
#endif

    //
    // We want a derived type that corresponds to the number of particles in
    // each cell. We only intend to use it in plotfiles for debugging purposes.
    // We'll just use the DERNULL since don't do anything in fortran for now.
    // We'll actually set the values in `writePlotFile()`.
    //
    derive_lst.add("particle_count", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(CA_DERNULL, ca_dernull), the_same_box);
    derive_lst.addComponent("particle_count", desc_lst, State_Type, Density, 1);

    derive_lst.add("particle_mass_density", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(CA_DERNULL, ca_dernull), grow_box_by_one);
    derive_lst.addComponent("particle_mass_density", desc_lst, State_Type,
                            Density, 1);

    derive_lst.add("total_particle_count", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(CA_DERNULL, ca_dernull), the_same_box);
    derive_lst.addComponent("total_particle_count", desc_lst, State_Type,
                            Density, 1);

    derive_lst.add("total_density", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(CA_DERNULL, ca_dernull), grow_box_by_one);
    derive_lst.addComponent("total_density", desc_lst, State_Type,
                            Density, 1);

    for (int i = 0; i < NumSpec; i++)
    {
        derive_lst.add(spec_names[i], IndexType::TheCellType(), 1,
                       BL_FORT_PROC_CALL(CA_DERSPEC, ca_derspec), the_same_box);
        derive_lst.addComponent(spec_names[i], desc_lst, State_Type, Density, 1);
        derive_lst.addComponent(spec_names[i], desc_lst, State_Type,
                                FirstSpec + i, 1);
    }

    for (int i = 0; i < NumAux; i++)
    {
        derive_lst.add(aux_names[i], IndexType::TheCellType(), 1,
                       BL_FORT_PROC_CALL(CA_DERSPEC, ca_derspec), the_same_box);
        derive_lst.addComponent(aux_names[i], desc_lst, State_Type, Density, 1);
        derive_lst.addComponent(aux_names[i], desc_lst, State_Type, FirstAux+i, 1);
    }
}

void
Nyx::no_hydro_setup()
{
    int dm = BL_SPACEDIM;

    //
    // Set number of state variables and pointers to components
    //
    int use_sgs = 0;
    int use_ssfm = 0;

    Density = 0;
    NUM_STATE = 1;

    // Define NUM_GROW from the f90 module.
    BL_FORT_PROC_CALL(GET_METHOD_PARAMS, get_method_params)(&NUM_GROW);

    BL_FORT_PROC_CALL(SET_METHOD_PARAMS, set_method_params)
        (dm, NumAdv, do_hydro, ppm_type, use_colglaz, 
         gamma, normalize_species ,use_sgs, use_ssfm);

    int coord_type = Geometry::Coord();
    BL_FORT_PROC_CALL(SET_PROBLEM_PARAMS, set_problem_params)
         (dm, phys_bc.lo(), phys_bc.hi(), Outflow, Symmetry, coord_type);

    Interpolater* interp = &cell_cons_interp;

    // Note that the default is state_data_extrap = false,
    // store_in_checkpoint = true.  We only need to put these in
    // explicitly if we want to do something different,
    // like not store the state data in a checkpoint directory
    bool state_data_extrap = false;
    bool store_in_checkpoint;

    store_in_checkpoint = true;
    desc_lst.addDescriptor(State_Type, IndexType::TheCellType(),
                           StateDescriptor::Point, 1, NUM_STATE, interp,
                           state_data_extrap, store_in_checkpoint);

#ifdef GRAVITY
    store_in_checkpoint = true;
    desc_lst.addDescriptor(PhiGrav_Type, IndexType::TheCellType(),
                           StateDescriptor::Point, 1, 1,
                           &cell_cons_interp, state_data_extrap,
                           store_in_checkpoint);
    store_in_checkpoint = false;
    desc_lst.addDescriptor(Gravity_Type, IndexType::TheCellType(),
                           StateDescriptor::Point, 1, BL_SPACEDIM,
                           &cell_cons_interp, state_data_extrap,
                           store_in_checkpoint);
#endif

    BCRec bc;
    set_scalar_bc(bc, phys_bc); 
    desc_lst.setComponent(State_Type, Density, "density", bc,
                          BndryFunc(BL_FORT_PROC_CALL(CA_DENFILL, ca_denfill)));

#ifdef GRAVITY
    if (do_grav)
    {
        set_scalar_bc(bc, phys_bc);
        desc_lst.setComponent(PhiGrav_Type, 0, "phi_grav", bc,
                              BndryFunc(BL_FORT_PROC_CALL(CA_PHIFILL,
                                                          ca_phifill)));
        set_x_vel_bc(bc, phys_bc);
        desc_lst.setComponent(Gravity_Type, 0, "grav_x", bc,
                              BndryFunc(BL_FORT_PROC_CALL(CA_GRAVXFILL,
                                                          ca_gravxfill)));
       set_y_vel_bc(bc, phys_bc);
       desc_lst.setComponent(Gravity_Type, 1, "grav_y", bc,
                             BndryFunc(BL_FORT_PROC_CALL(CA_GRAVYFILL,
                                                         ca_gravyfill)));
       set_z_vel_bc(bc, phys_bc);
       desc_lst.setComponent(Gravity_Type, 2, "grav_z", bc,
                             BndryFunc(BL_FORT_PROC_CALL(CA_GRAVZFILL,
                                                         ca_gravzfill)));
    }
#endif

#ifdef GRAVITY
    derive_lst.add("maggrav", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(CA_DERMAGGRAV, ca_dermaggrav),
                   the_same_box);
    derive_lst.addComponent("maggrav", desc_lst, Gravity_Type, 0, BL_SPACEDIM);
#endif

    //
    // We want a derived type that corresponds to the number of particles in
    // each cell. We only intend to use it in plotfiles for debugging purposes.
    // We'll just use the DERNULL since don't do anything in fortran for now.
    // We'll actually set the values in `writePlotFile()`.
    //
    derive_lst.add("particle_count", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(CA_DERNULL, ca_dernull), the_same_box);
    derive_lst.addComponent("particle_count", desc_lst, State_Type, Density, 1);

    derive_lst.add("particle_mass_density", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(CA_DERNULL, ca_dernull), grow_box_by_one);
    derive_lst.addComponent("particle_mass_density", desc_lst, State_Type,
                            Density, 1);

    derive_lst.add("total_particle_count", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(CA_DERNULL, ca_dernull), the_same_box);
    derive_lst.addComponent("total_particle_count", desc_lst, State_Type,
                            Density, 1);

    derive_lst.add("total_density", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(CA_DERNULL, ca_dernull), grow_box_by_one);
    derive_lst.addComponent("total_density", desc_lst, State_Type,
                            Density, 1);
}
