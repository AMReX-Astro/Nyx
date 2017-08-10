#include <iomanip>
#include "Nyx.H"

using namespace amrex;

void
Nyx::write_info ()
{
    int ndatalogs = parent->NumDataLogs();

    if (ndatalogs > 0)
    {
#ifndef NO_HYDRO
        MultiFab& D_new = get_new_data(DiagEOS_Type);
	Real      max_t = 0;

        Real rho_T_avg=0.0, T_avg=0.0, T_meanrho=0.0;
	if (do_hydro)
        {
            compute_new_temp();
            max_t = D_new.norm0(Temp_comp);
            compute_rho_temp(rho_T_avg, T_avg, T_meanrho);
	}
#endif

#ifdef NO_HYDRO
        Real time  = state[PhiGrav_Type].curTime();
#else
        Real time  = state[State_Type].curTime();
#endif
        Real dt    = parent->dtLevel(0);
        int  nstep = parent->levelSteps(0);

        if (ParallelDescriptor::IOProcessor())
        {
            std::ostream& data_loga = parent->DataLog(0);

            if (time == 0.0)
            {
                data_loga << std::setw( 8) <<  "   nstep";
                data_loga << std::setw(14) <<  "       time    ";
                data_loga << std::setw(14) <<  "        dt     ";
                data_loga << std::setw(14) <<  "      redshift ";
                data_loga << std::setw(14) <<  "       a       ";
#ifndef NO_HYDRO
                if (do_hydro == 1)
                {
                   data_loga << std::setw(14) <<  "  max temp     ";
                   data_loga << std::setw(14) <<  "rho-wgted temp ";
                   data_loga << std::setw(14) <<  " V-wgted temp  ";
                   data_loga << std::setw(14) <<  " T @ <rho>     ";
                }
#endif
                data_loga << '\n';

                Real old_z = (1. / old_a) - 1.;
                data_loga << std::setw( 8) <<  nstep;
                data_loga << std::setw(14) <<  std::setprecision(6) <<  time;
                data_loga << std::setw(14) <<  std::setprecision(6) <<    dt;
                data_loga << std::setw(14) <<  std::setprecision(6) << old_z;
                data_loga << std::setw(14) <<  std::setprecision(6) << old_a;
#ifndef NO_HYDRO
                if (do_hydro == 1)
                {
                   data_loga << std::setw(14) <<  std::setprecision(6) << max_t;
                   data_loga << std::setw(14) <<  std::setprecision(6) << rho_T_avg;
                   data_loga << std::setw(14) <<  std::setprecision(6) << T_avg;
                   data_loga << std::setw(14) <<  std::setprecision(6) << T_meanrho;
                }
#endif
                data_loga << '\n';
            }
            else
            {
                const Real new_z = (1. / new_a) - 1.;
                data_loga << std::setw( 8) <<  nstep;
                data_loga << std::setw(14) <<  std::setprecision(6) <<  time;
                data_loga << std::setw(14) <<  std::setprecision(6) <<    dt;
                data_loga << std::setw(14) <<  std::setprecision(6) << new_z;
                data_loga << std::setw(14) <<  std::setprecision(6) << new_a;
#ifndef NO_HYDRO
                if (do_hydro == 1)
                {
                   data_loga << std::setw(14) <<  std::setprecision(6) << max_t;
                   data_loga << std::setw(14) <<  std::setprecision(6) << rho_T_avg;
                   data_loga << std::setw(14) <<  std::setprecision(6) << T_avg;
                   data_loga << std::setw(14) <<  std::setprecision(6) << T_meanrho;
                }
#endif
                data_loga << std::endl;
            }
        }
    }
}
