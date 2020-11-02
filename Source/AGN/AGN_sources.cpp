#include <Nyx.H>

using namespace amrex;

#ifdef AGN
void
Nyx::get_old_source (Real      old_time,
                     Real      dt,
                     MultiFab& ext_src)
{
    BL_PROFILE("Nyx::get_old_source()");
    const Real* dx      = geom.CellSize();
    const Real* prob_lo = geom.ProbLo();
    const Real a        = get_comoving_a(old_time);
    const Real z        = 1 / a - 1;

    MultiFab& S_old = get_old_data(State_Type);
    MultiFab& D_old = get_old_data(DiagEOS_Type);

    ext_src.setVal(0.);

#if 0
    // Find the current particle locations
    Vector<Real> part_locs_and_mass;
    Nyx::theAPC()->GetParticleLocationsAndMass(part_locs_and_mass);
 
    Vector<Real> part_data;
    Nyx::theAPC()->GetParticleData(part_data);

    for (FillPatchIterator 
         Old_fpi (*this, S_old, 4, old_time, State_Type  , Density, S_old.nComp()),
         Old_dfpi(*this, D_old, 4, old_time, DiagEOS_Type, 0      , D_old.nComp());
         Old_fpi.isValid();
         ++Old_fpi)
    {
        const Box& bx = grids[Old_fpi.index()];
        ext_src
            (bx.loVect(), bx.hiVect(), 
             BL_TO_FORTRAN(Old_fpi()), BL_TO_FORTRAN(Old_fpi()), 
             BL_TO_FORTRAN(Old_dfpi()), BL_TO_FORTRAN(Old_dfpi()),
             BL_TO_FORTRAN(ext_src[Old_fpi]),
             part_locs_and_mass.dataPtr(), part_data.dataPtr(),
             prob_lo, dx, &old_time, &z, &dt);

        // The formulae in subroutine ctoprim assume that the source term for density is zero
        // Here we abort if it is non-zero.
        if (ext_src[Old_fpi].norm(0,Density,1) != 0)
        {
            std::cout << "The source terms for density are non-zero" << std::endl;
            amrex::Error();
        }
    }

    ext_src.EnforcePeriodicity(0, NUM_STATE, geom.periodicity());
#endif
}

void
Nyx::get_new_source (Real      old_time,
                     Real      new_time,
                     Real      dt,
                     MultiFab& ext_src)
{
    BL_PROFILE("Nyx::get_new_source()");
    const Real* dx      = geom.CellSize();
    const Real* prob_lo = geom.ProbLo();
    const Real a        = get_comoving_a(new_time);
    const Real z        = 1 / a - 1;

    MultiFab& S_old = get_old_data(State_Type);
    MultiFab& D_old = get_old_data(DiagEOS_Type);

    ext_src.setVal(0.);

#if 0
    // Find the current particle locations
    Vector<Real> part_locs_and_mass;
    Nyx::theAPC()->GetParticleLocationsAndMass(part_locs_and_mass);
    std::cout << "AGN LOCS " << part_locs_and_mass[0] << " " << part_locs_and_mass[1] << " " 
                             << part_locs_and_mass[2] << " " << part_locs_and_mass[3] << std::endl;
 
    Vector<Real> part_data;
    Nyx::theAPC()->GetParticleData(part_data);
    std::cout << "AGN DATA(V) " << part_data[0] << " " << part_data[1] << " " << part_data[2] << std::endl;
    std::cout << "AGN DATA(A) " << part_data[3] << " " << part_data[4] << " " << part_data[5] << std::endl;

    for (FillPatchIterator Old_fpi( *this, S_old, 4, old_time, State_Type  , Density, S_old.nComp()),
                           New_fpi( *this, S_old, 4, new_time, State_Type  , Density, S_old.nComp()),
                           Old_dfpi(*this, D_old, 4, old_time, DiagEOS_Type, 0      , D_old.nComp()),
                           New_dfpi(*this, D_old, 4, new_time, DiagEOS_Type, 0      , D_old.nComp());
         Old_fpi.isValid() && New_fpi.isValid() && Old_dfpi.isValid() && New_dfpi.isValid();
         ++Old_fpi, ++New_fpi, ++Old_dfpi, ++New_dfpi)
    {
        const Box& bx = grids[Old_fpi.index()];
        ext_src
            (bx.loVect(), bx.hiVect(), 
             BL_TO_FORTRAN(Old_fpi()), BL_TO_FORTRAN(New_fpi()), 
             BL_TO_FORTRAN(Old_dfpi()), BL_TO_FORTRAN(New_dfpi()), 
             BL_TO_FORTRAN(ext_src[Old_fpi]),
             part_locs_and_mass.dataPtr(), part_data.dataPtr(),
             prob_lo, dx, &new_time, &z, &dt);
    }

    ext_src.EnforcePeriodicity(0, NUM_STATE, geom.periodicity());
#endif
}
#endif
