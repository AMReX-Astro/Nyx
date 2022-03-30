#include <iomanip>
#include <Nyx.H>
#include <constants_cosmo.H>

using namespace amrex;

Real Nyx::initial_z             = -1.0;
Real Nyx::final_a               = -1.0;
Real Nyx::final_z               = -1.0;
Real Nyx::relative_max_change_a =  0.01;
Real Nyx::absolute_max_change_a = -1.0;
Real Nyx::dt_binpow             = -1.0;
Real Nyx::initial_time          = -1.0;
Real Nyx::final_time            = -1.0;

void
Nyx::read_comoving_params ()
{
    ParmParse pp("nyx");

    pp.get("comoving_OmB"   , comoving_OmB);
    pp.get("comoving_OmM"   , comoving_OmM);
    pp.get("comoving_h"     , comoving_h);
    pp.query("comoving_OmR" , comoving_OmR);

    pp.query("initial_z", initial_z);
    pp.query("final_a",   final_a);
    pp.query("final_z",   final_z);

    if (final_z >= 0)
    {
        if (final_a > 0)
        {
            std::cerr << "ERROR::dont specify both final_a and final_z\n";
            amrex::Error();
        }
        else
        {
            final_a = 1 / (1 + final_z);
        }
    }

    pp.query("relative_max_change_a", relative_max_change_a);
    pp.query("absolute_max_change_a", absolute_max_change_a);
    pp.query("dt_binpow",             dt_binpow);


    // for shrinking box tests, initial_z < 0 is ok
    if (initial_z < 0)
    {
        std::cerr << "ERROR::Need to specify non-negative initial redshift \n";
        amrex::Error();
    }

    if (comoving_h > 0)
    {
       // save start/end times, for reporting purposes
       Real a0 = 0.0, a1 = 1.0/(1.0+initial_z);
       integrate_time_given_a(a0, a1, initial_time);

       if (final_z >= 0) {
          a1 = 1.0/(1.0+final_z);
       } else {
          a1 = final_a;
       }
       integrate_time_given_a(a0, a1, final_time);
    } else {
       // These are just defaults so the values are defined
       initial_time = 0.0;
         final_time = 0.0;
    }
}

void
Nyx::integrate_comoving_a_to_a(const Real old_a_local, const Real a_value, Real& dt)
{
    Real H_0, OmL;
    Real Delta_t;
    Real start_a, end_a, start_slope, end_slope;
    int j, nsteps;

    if (comoving_h == 0.0)
        amrex::Abort("fort_integrate_comoving_a_to_z: Shouldn't be setting plot_z_values if not evolving a");

    H_0 = comoving_h * Hubble_const;
    OmL = 1.e0 - comoving_OmM - comoving_OmR;

    // Use lots of steps if we want to nail the z_value
    //        nsteps = 1024

    // Use enough steps if we want to be close to the a_value
    nsteps = 2048;

    // We integrate a, but stop when a = a_value (or close enough)
    Delta_t = dt/nsteps;
    end_a = old_a_local;
    for(j = 1; j <= nsteps; j++)
    {
        // This uses RK2 to integrate the ODE:
        //   da / dt = H_0 * sqrt(OmM/a + OmR/a^2 + OmL*a^2)
        start_a = end_a;

        // Compute the slope at the old time
        if (comoving_type > 0)
            start_slope = H_0*std::sqrt(comoving_OmM/start_a + comoving_OmR/(start_a*start_a) + OmL*start_a*start_a);
        else
            start_slope = comoving_h;

        // Compute a provisional value of ln(a) at the new time
        end_a = start_a + start_slope * Delta_t;

        // Compute the slope at the new time
        if (comoving_type > 0)
            end_slope = H_0*std::sqrt(comoving_OmM/end_a + comoving_OmR/(end_a*end_a) + OmL*end_a*end_a);
        else
            end_slope = comoving_h;

        // Now recompute a at the new time using the average of the two slopes
        end_a = start_a + 0.5e0 * (start_slope + end_slope) * Delta_t;

        // We have crossed from a too small to a too big in this step
        if ( (end_a - a_value) * (start_a - a_value) < 0)
        {
            dt = ( (  end_a - a_value) * double(j  ) +
                   (a_value - start_a) * double(j+1) ) / (end_a - start_a) * Delta_t;
            break;
        }

     }
}

void
Nyx::integrate_comoving_a_to_z(const Real old_a_local, const Real z_value, Real& dt)
{
    Real H_0, OmL;
    Real Delta_t;
    Real start_a, end_a, start_slope, end_slope;
        Real a_value;
    int j, nsteps;

    if (comoving_h == 0.0)
        amrex::Abort("fort_integrate_comoving_a_to_z: Shouldn't be setting plot_z_values if not evolving a");

    H_0 = comoving_h * Hubble_const;
    OmL = 1.e0 - comoving_OmM - comoving_OmR;

    // Translate the target "z" into a target "a"
    a_value = 1.e0 / (1.e0 + z_value);

    // Use lots of steps if we want to nail the z_value
    nsteps = 1024;

    // We integrate a, but stop when a = a_value (or close enough)
    Delta_t = dt/nsteps;
    end_a = old_a_local;
    for(j = 1; j <= nsteps; j++)
    {
        // This uses RK2 to integrate the ODE:
        //   da / dt = H_0 * sqrt(OmM/a + OmR/a^2 + OmL*a^2)
        start_a = end_a;

        // Compute the slope at the old time
        if (comoving_type > 0)
            start_slope = H_0*std::sqrt(comoving_OmM/start_a + comoving_OmR/(start_a*start_a) + OmL*start_a*start_a);
        else
            start_slope = comoving_h;

        // Compute a provisional value of ln(a) at the new time
        end_a = start_a + start_slope * Delta_t;

        // Compute the slope at the new time
        if (comoving_type > 0)
            end_slope = H_0*std::sqrt(comoving_OmM/end_a + comoving_OmR/(end_a*end_a) + OmL*end_a*end_a);
        else
            end_slope = comoving_h;

        // Now recompute a at the new time using the average of the two slopes
        end_a = start_a + 0.5e0 * (start_slope + end_slope) * Delta_t;

        // We have crossed from a too small to a too big in this step
        if ( (end_a - a_value) * (start_a - a_value) < 0)
        {
            dt = ( (  end_a - a_value) * double(j  ) +
                   (a_value - start_a) * double(j+1) ) / (end_a - start_a) * Delta_t;
            break;
        }

     }
}

void
Nyx::est_maxdt_comoving_a(const Real old_a_local, Real & dt)
{
    Real H_0, OmL;
    Real max_dt;

    OmL = 1.e0 - comoving_OmM - comoving_OmR;

    // This subroutine computes dt based on not changing a by more than 5%
    // if we use forward Euler integration
    //   d(ln(a)) / dt = H_0 * sqrt(OmM/a^3 + OmL)

    H_0 = comoving_h * Hubble_const;

    if (H_0 != 0.0)
    {
        if (comoving_type > 0)
            max_dt = (0.05) / H_0 / std::sqrt( comoving_OmM/(old_a_local * old_a_local *old_a_local) 
                                              +comoving_OmR/(old_a_local * old_a_local *old_a_local * old_a_local) + OmL);
        else
            max_dt = (0.05) / std::abs(comoving_h);
        dt = amrex::min(dt,max_dt);
    }
}

// This might only work for t=0=> a=.00625, although constant canceled
void
Nyx::est_lindt_comoving_a(const Real old_a_local, const Real new_a_local, Real& dt)
{
    Real H_0;
    Real lin_dt;

    // This subroutine computes dt based on not changing a by more than 5%
    // if we use forward Euler integration
    // OmL = 1.e0 - comoving_OmM - comoving_OmR;
    //   d(ln(a)) / dt = H_0 * sqrt(OmM/a^3 + OmR/a^4 + OmL)

    H_0 = comoving_h * Hubble_const;

    // Could definately be optimized better
    if (H_0 != 0.0)
    {
        lin_dt= ( pow((new_a_local/(pow(.75,(2_rt/3_rt)))),(.75))
                 -pow((old_a_local/(pow(.75,(2_rt/3_rt)))),(.75))  ) / H_0;
        dt=lin_dt;
    }
}

void
Nyx::estdt_comoving_a(const Real old_a_local, Real& new_a_local, 
                      Real& dt, const Real change_allowed, 
                      const Real fixed_da, const Real final_a_in, int& dt_modified)
{
    if (comoving_h != 0.0e0)
    {
        if( fixed_da <= 0.0e0)
        {
            // First call this to make sure dt that we send to integration routine isnt outrageous
            est_maxdt_comoving_a(old_a_local,dt);

            // Initial call to see if existing dt will work
            integrate_comoving_a(old_a_local,new_a_local,dt);

            // Make sure a isn't growing too fast
            enforce_percent_change(old_a_local,new_a_local,dt,change_allowed);
        }
        else
        {
            // First call this to make sure dt that we send to integration routine isnt outrageous
            new_a_local = (old_a_local +  fixed_da);
            est_lindt_comoving_a(old_a_local,new_a_local,dt);
            est_maxdt_comoving_a(old_a_local,dt);

            // Then integrate old_a_local to a_value using dt as a guess for the maximum dt
            // output dt is based on a fraction of the input dt
            integrate_comoving_a_to_a(old_a_local,new_a_local,dt);
        }
        // Make sure we don't go past final_a_in (if final_a_in is set)
        if (final_a_in > 0.0e0)
            enforce_final_a(old_a_local,new_a_local,dt,final_a_in);

        dt_modified = 1;
    }
    else
    {
        // dt is unchanged by this call

        dt_modified = 0;
    }
}

void
Nyx::enforce_percent_change(const Real old_a_local, Real& new_a_local, 
                            Real& dt,const Real change_allowed)
{

    int i;
    Real factor = ( (new_a_local - old_a_local) / old_a_local ) / change_allowed;

    // Only go into this process if percent change exceeds change_allowed

    if(factor > 1.0)
    {
        for(i = 1; i <= 100; i++)
        {
            factor = ( (new_a_local - old_a_local) / old_a_local ) / change_allowed;

            // Note: 0.99 is just a fudge factor so we don't get bogged down.
            if(factor > 1.0)
            {
                dt = (1.0 / factor) * dt * 0.99;
                integrate_comoving_a(old_a_local,new_a_local,dt);
            }
            else if (i < 100)
            {
                integrate_comoving_a(old_a_local,new_a_local,dt);
                // We're done
                return;
            }
            else
                amrex::Abort("Too many iterations in enforce_percent_change");

        }
    }
    else
        return;
}

void
Nyx::enforce_final_a(const Real old_a_local, Real& new_a_local,
                     Real& dt, const Real final_a_in)
{
    int i;
    Real factor;
    const Real eps = 1.e-10;

    if (old_a_local > final_a_in)
        amrex::Abort("Oops -- old_a > final_a_in");

    // Only go into this process if new_a is past final_a_in
    if (new_a_local > final_a_in)
    {
        for(i = 1; i <= 100; i++)
        {
            if ( (new_a_local > (final_a_in+eps)) || (new_a_local < final_a_in) )
            {
                factor = (final_a_in - old_a_local) / (new_a_local - old_a_local);
                dt = dt * factor;
                integrate_comoving_a(old_a_local,new_a_local,dt);
            }
            else if (i<100)
                return;
            else
                amrex::Abort("Too many iterations in enforce_final_a");
        }
    }
    else
        // We don't need to do anything
        return;
}


void
Nyx::comoving_est_time_step (Real& cur_time, Real& estdt)
{
    Real change_allowed = relative_max_change_a;
    Real fixed_da = absolute_max_change_a;
    Real dt             = estdt;
    Real new_dummy_a;
    int  dt_modified;

    if ( std::abs(cur_time - new_a_time) <= 1.e-12 * cur_time)
    {

        // Initial guess -- note that we send in "new_a" because we haven't yet swapped
        // "old_a" and "new_a" -- we can't do that until after we compute dt and then
        // integrate a forward.
        estdt_comoving_a
          (new_a, new_dummy_a, dt, change_allowed, fixed_da, final_a, dt_modified);

    if(dt_binpow >= 0)
      {
        if(estdt>=dt)
          estdt=dt;
        else if(estdt>.5*dt)
          {
            estdt=.5*dt;
            //      std::cout << "Lavel = 1" <<std::endl;
          }
        else if(estdt>.25*dt)
          {
            estdt=.25*dt;
            //      std::cout << "Lavel = 2" <<std::endl;
          }
        else if(estdt>.125*dt)
          {
            estdt=.125*dt;
            //      std::cout << "Lavel = 3" <<std::endl;
          }
        else if(estdt>.0625*dt)
          {
            estdt=.0625*dt;
            //      std::cout << "Lavel = 4" <<std::endl;
          }
        else
          {
            //dta*(2**(-1*np.ceil( np.log2(dta/dth))))
            estdt = dt*(pow(2,(-std::ceil( std::log2(dt/estdt)))));
            //      std::cout << "Lavel > 4" <<std::endl;
          }
        integrate_comoving_a(new_a,new_dummy_a,estdt);
      }
    else
      {
        estdt=std::min(estdt,dt);
      }
          

        if (verbose && (dt_modified == 1) && ParallelDescriptor::IOProcessor())
        {
            std::cout << "...estdt after call to comoving   : "
                      << dt
                      << "\n...change in a is "
                      << (new_dummy_a - new_a) / new_a * 100.0
                      << " percent\n";
        }
    } 
    // This only happens with the second call at time = 0 where we want to re-advance
    //      from old_a_time time, not start at new_a_time
    else if ( std::abs(cur_time - old_a_time) <= 1.e-12 * cur_time)
    {

        // Initial guess -- note that we send in "new_a" because we haven't yet swapped
        // "old_a" and "new_a" -- we can't do that until after we compute dt and then
        // integrate a forward.
        estdt_comoving_a
          (old_a, new_dummy_a, dt, change_allowed, fixed_da, final_a, dt_modified);
    if(dt_binpow >= 0)
      {
        if(estdt>=dt)
          estdt=dt;
        else if(estdt>.5*dt)
          {
            estdt=.5*dt;
            //      std::cout << "Lavel = 1" <<std::endl;
          }
        else if(estdt>.25*dt)
          {
            estdt=.25*dt;
            //      std::cout << "Lavel = 2" <<std::endl;
          }
        else if(estdt>.125*dt)
          {
            estdt=.125*dt;
            //      std::cout << "Lavel = 3" <<std::endl;
          }
        else if(estdt>.0625*dt)
          {
            estdt=.0625*dt;
            //      std::cout << "Lavel = 4" <<std::endl;
          }
        else
          {
            //dta*(2**(-1*np.ceil( np.log2(dta/dth))))
            estdt = dt*(pow(2,(-std::ceil( std::log2(dt/estdt)))));
            //      std::cout << "Lavel > 4" <<std::endl;
          }
        integrate_comoving_a(old_a,new_dummy_a,estdt);
      }
    else
      {
        estdt=std::min(estdt,dt);
      }

        if (verbose && (dt_modified == 1) && ParallelDescriptor::IOProcessor())
        {
            std::cout << "...advancing from old_a_time rather than new_a_time! " << std::endl;
            std::cout << "...estdt after call to comoving   : "
                      << dt
                      << "\n...change in a is "
                      << (new_dummy_a - old_a) / old_a * 100.0
                      << " percent\n";
        }

    } 
    else if ( std::abs(cur_time - old_a_time) <= 1.e-12 * cur_time)
    {
       std::cout << "comoving_est_time_step: DONT KNOW WHAT TIME IT IS " << cur_time << std::endl;
       exit(0);
    } 

    return;
}

Real
Nyx::get_comoving_a (Real time)
{
    const Real eps         = 0.0001 * (new_a_time - old_a_time);

    // Test on whether time == old_a_time == new_a_time -- for example after restart
    //   before a has been integrated for the new time step.
    if ( ( std::abs(time - old_a_time) <= 1.e-12*old_a_time ) &&
         ( std::abs(time - new_a_time) <= 1.e-12*new_a_time ) &&
         (  old_a == new_a ) )
    {
        return old_a;
    }
    else if (time > old_a_time - eps && time < old_a_time + eps)
    {
        return old_a;
    }
    else if (time > new_a_time - eps && time < new_a_time + eps)
    {
        return new_a;
    }
    else if (time > old_a_time && time < new_a_time)
    {
        Real frac = (time - old_a_time) / (new_a_time - old_a_time);
        Real    a = frac*new_a + (1.0-frac)*old_a;
        return a;
    } 
    else
    {
        if (verbose && ParallelDescriptor::IOProcessor())
        {
            std::cout << "Invalid:get comoving_a at " << time << std::endl; 
            std::cout << "Old / new a_time  " << old_a_time << " " << new_a_time << std::endl;
            std::cout << "Old / new   a     " << old_a      << " " << new_a    << std::endl;
        }
        amrex::Error("get_comoving_a: invalid time");
        return 0;
    }
}

void
Nyx::integrate_comoving_a (Real time,Real dt)
{
    if (level > 0) 
        return;

    bool first;

    if ( std::abs(time-new_a_time) <= (1.e-10 * time) )
    {
       first = true;
    } else {
       first = false;
    } 

    if (first) 
    {

        // Update a
        old_a      = new_a;
        integrate_comoving_a(old_a, new_a, dt);

        // Update the times
        old_a_time = new_a_time;
        new_a_time = old_a_time + dt;
 
        if (verbose && ParallelDescriptor::IOProcessor())
        {
            std::cout << "Integrating a from time " << time << " by dt = " << dt << '\n';
            std::cout << "Old / new A time      " << old_a_time << " " << new_a_time << std::endl;
            std::cout << "Old / new A           " << old_a      << " " << new_a      << std::endl;
            std::cout << "Old / new z           " << 1./old_a-1.<< " " << 1./new_a-1. << std::endl;
        }
    }
    else if (std::abs(time-old_a_time) <= 1.e-10 * time) 
    {
        // Leave old_a and old_a_time alone -- we have already swapped them
        integrate_comoving_a(old_a, new_a, dt);

        // Update the new time only
        new_a_time = old_a_time + dt;

        if (verbose && ParallelDescriptor::IOProcessor())
        {
            std::cout << "Re-integrating a from time " << time << " by dt = " << dt << '\n';
            std::cout << "Old / new A time         " << old_a_time << " " << new_a_time << std::endl;
            std::cout << "Old / new A              " << old_a      << " " << new_a      << std::endl;
            std::cout << "Old / new z              " << 1./old_a-1.<< " " << 1./new_a-1. << std::endl;
        }
    }
    else 
    {
            std::cout << "Time passed to integrate_comoving_a " << time << std::endl;
            std::cout << "Old / new A time                    " << old_a_time << " " << new_a_time << std::endl;
            std::cout << "Old / new A                         " << old_a      << " " << new_a      << std::endl;
            std::cerr << "ERROR::dont know what to do in integrate_comoving_a" << std::endl;
            amrex::Error();
    }
}

Real invEz(Real H0, Real Om, Real Or, Real a)
{
        //      invEz = 
        return 1.0e0 / ( H0*std::sqrt(Om/a + Or/(a*a) + (1.0e0-Om-Or)*a*a) );
}

void
Nyx::integrate_time_given_a(const Real a0, const Real a1, Real& dt)
{
    const Real small_dt_fac = 1.0e-6;
    Real H0, OmM, OmR, prev_soln, h;
    int iter, n, j;

    H0  = comoving_h*Hubble_const;
    OmM = comoving_OmM;
    OmR = comoving_OmR;

    prev_soln = -1.0e0;
    // trapezoidal integration
    for(iter = 1; iter<= 20; iter++)
    {
        n  = static_cast<int>(std::round(std::pow(2,iter)));

        h = (a1-a0)/(n-1);

        if (a0 < 1.0e-10) // prevent division by zero in invEz
            dt = 0.5*invEz(H0, OmM, OmR, a1);
        else
            dt = 0.5*(invEz(H0, OmM, OmR, a0) + invEz(H0, OmM, OmR, a1));

        for(j = 1; j <= n-2; j++)
        {
            dt = dt + invEz(H0, OmM, OmR, a0+j*h);
        }
        dt = dt*h;

        if (iter > 4)
            if (std::abs(dt-prev_soln) < small_dt_fac*std::abs(prev_soln))
                return;

        prev_soln = dt;
    }
}

void
Nyx::integrate_comoving_a (const Real old_a_local, Real& new_a_local, const Real dt)
{

    const Real small_a_fac = 1.0e-8;
    Real H_0, OmL, Delta_t, prev_soln;
    Real start_a, end_a, start_slope, end_slope;
    int iter, j, nsteps;

    if (comoving_h == 0.0)
    {
        new_a_local = old_a_local;
        return;
    }

    H_0 = comoving_h * Hubble_const;
    OmL = 1.e0 - comoving_OmM - comoving_OmR;

    prev_soln = 2.0e0; // 0<a<1 so a=2 will do as "wrong" solution
    for(iter = 1; iter<=11; iter++)//  ! max allowed iterations
    {
        nsteps  = static_cast<int>(std::round(std::pow(2,iter-1)));

        Delta_t = dt/nsteps;
        end_a = old_a_local;

        for(j = 1; j <= nsteps; j++)
        {
            // This uses RK2 to integrate the ODE:
            //   da / dt = H_0 * sqrt(OmM/a + OmR/a^2 + OmL*a^2)
            start_a = end_a;

            // Compute the slope at the old time
            if (comoving_type > 0)
                start_slope = H_0*std::sqrt(comoving_OmM/start_a + comoving_OmR/(start_a*start_a) + OmL*start_a*start_a);
            else
                start_slope = comoving_h;

            // Compute a provisional value of ln(a) at the new time
            end_a = start_a + start_slope * Delta_t;

            // Compute the slope at the new time
            if (comoving_type > 0)
                end_slope = H_0*std::sqrt(comoving_OmM/end_a + comoving_OmR/(end_a*end_a) + OmL*end_a*end_a);
            else
                end_slope = comoving_h;

            // Now recompute a at the new time using the average of the two slopes
            end_a = start_a + 0.5e0 * (start_slope + end_slope) * Delta_t;

        }

        new_a_local  = end_a;

        if (std::abs(1.0e0-new_a_local/prev_soln) <= small_a_fac)
              return;
        prev_soln = new_a_local;

    }
}

void
Nyx::comoving_a_post_restart (const std::string& restart_file)
{
    if (level > 0)
        return;

    if (ParallelDescriptor::IOProcessor())
    {
        std::string FileName = restart_file + "/comoving_a";
        std::ifstream File;
        File.open(FileName.c_str(),std::ios::in);
        if (!File.good())
            amrex::FileOpenFailed(FileName);
        File >> old_a;
    }
    ParallelDescriptor::Bcast(&old_a, 1, ParallelDescriptor::IOProcessorNumber());

    new_a = old_a;

#ifdef NO_HYDRO
    old_a_time = state[PhiGrav_Type].curTime();
#else
    old_a_time = state[State_Type].curTime();
#endif
    new_a_time = old_a_time;
    if (ParallelDescriptor::IOProcessor())
    {
        std::cout << "...setting old_a_time to " << old_a_time << std::endl;
    }

}

void
Nyx::plot_z_est_time_step (Real& dt_0, bool& dt_changed)
{
    Real dt = dt_0;
    Real a_old, z_old, a_new, z_new;
    Real z_value;

    // This is where we are now
#ifdef NO_HYDRO
    Real cur_time = state[PhiGrav_Type].curTime();
#else
    Real cur_time = state[State_Type].curTime();
#endif
    a_old = get_comoving_a(cur_time);
    z_old = (1. / a_old) - 1.;

    // *****************************************************
    // First test whether we are within dt of a plot_z value
    // *****************************************************

    // This is where we would be if we use the current dt_0
    Real new_time = cur_time + dt_0;
    integrate_comoving_a (cur_time,dt_0);
    a_new = get_comoving_a(new_time);
    z_new = (1. / a_new) - 1.;

    // Find the relevant entry of the plot_z_values array
    bool found_one = false;
    for (int i = 0; i < plot_z_values.size(); i++)
    {
        // We have gone from before to after one of the specified values
        if ( (z_new - plot_z_values[i]) * (z_old - plot_z_values[i]) < 0 && !found_one)
        {
            z_value   = plot_z_values[i];
            found_one = true;
        }
    }

    // Now that we know that dt_0 is too big and makes us pass z_value,
    // we must figure out what value of dt < dt_0 makes us exactly reach z_value
    if (found_one)
    {
        integrate_comoving_a_to_z(old_a, z_value, dt);
                

        if (verbose && ParallelDescriptor::IOProcessor())
        {
            std::cout << " " << std::endl;
            std::cout << " ... modifying time step from " << dt_0 << " to " << dt << std::endl;
            std::cout << " ... in order to write a plotfile at z = " << z_value << std::endl;
            std::cout << " " << std::endl;
        }

        // We want to pass this value back out
        dt_0 = dt;
        dt_changed = true;
    } 
    else 
    { 

    // *****************************************************
    // If not within one dt, now test whether we are within 2*dt of a plot_z value
    // *****************************************************

    // This is where we would be if we advance by twice the current dt_0
    Real two_dt = 2.0*dt_0;

    integrate_comoving_a(a_old, a_new, two_dt);

    z_new = (1. / a_new) - 1.;

    // Find the relevant entry of the plot_z_values array
    found_one = false;
    for (int i = 0; i < plot_z_values.size(); i++)
    {
        // We have gone from before to after one of the specified values
        if ( (z_new - plot_z_values[i]) * (z_old - plot_z_values[i]) < 0 && !found_one)
        {
            z_value   = plot_z_values[i];
            found_one = true;
        }
    }

    // Now that we know that 2*dt_0 will make us pass z_value, we set the current dt
    // as half the interval to reach that z_value
    if (found_one)
    {
        integrate_comoving_a_to_z(old_a, z_value, two_dt);
                

        if (verbose && ParallelDescriptor::IOProcessor())
        {
            std::cout << " " << std::endl;
            std::cout << " ... modifying time step from " << dt_0 << " to " << 0.5 * two_dt << std::endl;
            std::cout << " ... in order to write a plotfile at z = " << z_value 
                      << " two steps from now " << std::endl;
            std::cout << " " << std::endl;
        }

        // We want to pass this value back out
        dt_0 = 0.5 * two_dt;
        dt_changed = true;
    }

    }
}

void
Nyx::analysis_z_est_time_step (Real& dt_0, bool& dt_changed)
{
    Real dt = dt_0;
    Real a_old, z_old, a_new, z_new, z_value;

    // This is where we are now
#ifdef NO_HYDRO
    Real cur_time = state[PhiGrav_Type].curTime();
#else
    Real cur_time = state[State_Type].curTime();
#endif
    a_old = get_comoving_a(cur_time);
    z_old = (1. / a_old) - 1.;

    // *****************************************************
    // First test whether we are within dt of a analysis_z value
    // *****************************************************

    // This is where we would be if we use the current dt_0
    Real new_time = cur_time + dt_0;
    integrate_comoving_a (cur_time,dt_0);
    a_new = get_comoving_a(new_time);
    z_new = (1. / a_new) - 1.;

    // Find the relevant entry of the analysis_z_values array
    bool found_one = false;
    for (int i = 0; i < analysis_z_values.size(); i++)
    {
        // We have gone from before to after one of the specified values
        if ( (z_new - analysis_z_values[i]) * (z_old - analysis_z_values[i]) < 0 && !found_one)
        {
            z_value   = analysis_z_values[i];
            found_one = true;
        }
    }

    // Now that we know that dt_0 is too big and makes us pass z_value,
    // we must figure out what value of dt < dt_0 makes us exactly reach z_value
    if (found_one)
    {
        integrate_comoving_a_to_z(old_a, z_value, dt);

        if (verbose && ParallelDescriptor::IOProcessor())
        {
            std::cout << " " << std::endl;
            std::cout << " ... modifying time step from " << dt_0 << " to " << dt << std::endl;
            std::cout << " ... in order to do analysis at z = " << z_value << std::endl;
            std::cout << " " << std::endl;
        }

        // We want to pass this value back out
        dt_0 = dt;
        dt_changed = true;
    } 
    else 
    { 

    // *****************************************************
    // If not within one dt, now test whether we are within 2*dt of a analysis_z value
    // *****************************************************

    // This is where we would be if we advance by twice the current dt_0
    Real two_dt = 2.0*dt_0;

    integrate_comoving_a(a_old, a_new, two_dt); 
    z_new = (1. / a_new) - 1.;

    // Find the relevant entry of the analysis_z_values array
    found_one = false;
    for (int i = 0; i < analysis_z_values.size(); i++)
    {
        // We have gone from before to after one of the specified values
        if ( (z_new - analysis_z_values[i]) * (z_old - analysis_z_values[i]) < 0 && !found_one)
        {
            z_value   = analysis_z_values[i];
            found_one = true;
        }
    }

    // Now that we know that 2*dt_0 will make us pass z_value, we set the current dt
    // as half the interval to reach that z_value
    if (found_one)
    {
        two_dt = 2.0*dt_0;
        integrate_comoving_a_to_z(old_a, z_value, two_dt);

        if (verbose && ParallelDescriptor::IOProcessor())
        {
            std::cout << " " << std::endl;
            std::cout << " ... modifying time step from " << dt_0 << " to " << 0.5 * two_dt << std::endl;
            std::cout << " ... in order to do analysis at z = " << z_value 
                      << " two steps from now " << std::endl;
            std::cout << " " << std::endl;
        }

        // We want to pass this value back out
        dt_0 = 0.5 * two_dt;
        dt_changed = true;
    }

    }
}
