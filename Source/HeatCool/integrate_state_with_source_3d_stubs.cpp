#include <Nyx.H>

int Nyx::integrate_state_struct
  (amrex::MultiFab &S_old,
   amrex::MultiFab &S_new,
   amrex::MultiFab &D_old,
   amrex::MultiFab &hydro_src,
   amrex::MultiFab &IR,
   amrex::MultiFab &reset_src,
   const Real& a, const amrex::Real& a_end,
   const Real& delta_time,
   const int sdc_iter)
{
   return -1;
}

#ifdef SAVE_REACT
int Nyx::integrate_state_struct_mfin
  (amrex::Array4<Real> const& state4,
   amrex::Array4<Real> const& diag_eos4,
   amrex::Array4<Real> const& state_n4,
   amrex::Array4<Real> const& hydro_src4,
   amrex::Array4<Real> const& reset_src4,
   amrex::Array4<Real> const& IR4,
   amrex::Array4<Real> const& react_in_arr,
   amrex::Array4<Real> const& react_out_arr,
   amrex::Array4<Real> const& react_out_work_arr,
   const Box& tbx,
   const Real& a, const amrex::Real& a_end, const Real& delta_time,
   long int& old_max_steps, long int& new_max_steps,
   const int sdc_iter)
{
        return -1;
}
#else
int Nyx::integrate_state_struct_mfin
  (amrex::Array4<Real> const& state4,
   amrex::Array4<Real> const& diag_eos4,
   amrex::Array4<Real> const& state_n4,
   amrex::Array4<Real> const& hydro_src4,
   amrex::Array4<Real> const& reset_src4,
   amrex::Array4<Real> const& IR4,
   const Box& tbx,
   const Real& a, const amrex::Real& a_end, const Real& delta_time,
   long int& old_max_steps, long int& new_max_steps,
   const int sdc_iter)
{
        return -1;
}
#endif
