
Units and Conventions
=====================

Nyx uses a cosmological system of units based on Mpc, M\ :math:`_\odot`, and km/s for length, mass, and velocity,
respectively.  Temperature is given in degrees Kelvin.  The unit of time is derived from velocity and length units.
All inputs and problem initialization should be specified consistently with these units,
and the outputs should be interpreted with these units in mind.
In the equation of state calls, there is an internal
conversion between cosmological and CGS units.

:numref:`table:units` shows 
some of the common symbols / names used throughout the code documentation and papers.

.. _table:units:

.. table:: Common quantities and units.

   +-----------------------+--------------------------------------------------------+-----------------------+
   | name                  | units                                                  | description           |
   +=======================+========================================================+=======================+
   | :math:`t`             | (Mpc/km) s                                             | time                  |
   +-----------------------+--------------------------------------------------------+-----------------------+
   | :math:`\rho`          | M\ :math:`_\odot` / Mpc\ :math:`^3`                    | mass density          |
   +-----------------------+--------------------------------------------------------+-----------------------+
   | :math:`\ub`           | km/s                                                   | velocity vector       |
   +-----------------------+--------------------------------------------------------+-----------------------+
   | :math:`p`             | M\ :math:`_\odot` (km/s)\ :math:`^2` / Mpc\ :math:`^3` | pressure              |
   |                       |                                                        |                       |
   +-----------------------+--------------------------------------------------------+-----------------------+
   | :math:`\gb`           | (km/s)2 / Mpc                                          | gravitational         |
   |                       |                                                        | acceleration          |
   +-----------------------+--------------------------------------------------------+-----------------------+
   | :math:`E`             | (km/s)\ :math:`^2`                                     | specific total energy |
   +-----------------------+--------------------------------------------------------+-----------------------+
   | :math:`e`             | (km/s)\ :math:`^2`                                     | specific internal     |
   |                       |                                                        | energy                |
   +-----------------------+--------------------------------------------------------+-----------------------+
   | :math:`T`             | :math:`K`                                              | temperature           |
   +-----------------------+--------------------------------------------------------+-----------------------+


In :numref:`table:constants` we list the values used for physical constants in cosmological units.
Note that :math:`\Omega_m`, :math:`\Omega_b`, :math:`\Omega_r`  and :math:`h` are set in the inputs file :ref:`Comoving: List of Parameters<list-of-parameters-13>`.
Full list of constants and conversion factors is set in Source/Driver/constants_cosmo.H.

.. _table:constants:
.. table::
	   Physical constant values

   +-------------------------------------+-----------------------------------------------------------+
   | Constant                            | Cosmological units                                        |       
   +=====================================+===========================================================+
   | Gravitational constant (:math:`G`)  | 4.3019425e-9 Mpc (km/s)\ :math:`^2` M\ :math:`_\odot^{-1}`|
   +-------------------------------------+-----------------------------------------------------------+
   | Avogadro’s number (:math:`n_A`)     | 1.1977558e57 M\ :math:`_\odot^{-1}`                       |
   +-------------------------------------+-----------------------------------------------------------+
   | Boltzmann’s constant (:math:`k_B`)  | 0.6941701e-59 M\ :math:`_\odot` (km/s)\ :math:`^2` / K    |
   +-------------------------------------+-----------------------------------------------------------+
   | Hubble constant (:math:`H`)         | 100 (km/s) / Mpc                                          |
   +-------------------------------------+-----------------------------------------------------------+

The only other place that dimensional numbers are used in the code is in the tracing and Riemann solve.
We set three *small* numbers which need to be consistent with the data specified
We set one *large* number which needs to be consistent with the data specified.
Each of these can be specified in the inputs file.

-  small_dens – small value for density

-  small_press – small value for pressure

-  small_temp – small value for temperature

-  large_temp – large value for temperature

These are the places that each is used in the code:

-  **small_dens**

   -  | **enforce_minimum_density** (called after subroutine consup) – there are two choices for this. In the flooring routine, 
      | **enforce_minimum_density_floor** – density is set to small_dens, (rho e) and (rho E) are computed from small_temp,
      | and momenta are set to zero.  In the conservative routine, **enforce_minimum_density_cons**, an iterative procedure 
      | is used to create diffusive fluxes that adjusts all the variables conservatively until density is greater than small_dens.

   -  | **tracexy / tracez / tracexy_ppm / tracez_ppm**:
      | qxp = max(qxp,small_dens)
      | qxm = max(qxm,small_dens)
      | and analogously for qyp/qym and qzp/qzm. This modifies the primitive density and pressure inside the tracing, not the underlying state variables

   -  **riemannus** – we set

      wsmall = small_dens \* csmall

      and then

      | wl = max(wsmall, sqrt(gaml \* pl \* rl))
      | wr = max(wsmall, sqrt(gamr \* pr \* rr))

      Also, we set

      ro = max(small_dens,ro)

      where ro = 0.5 \* (rl + rr) – this state is only chosen when ustar = 0, and

      rstar = max(small_dens,rstar)

      where rstar = ro + (pstar-po)/co:math:`^2`

-  **small_temp**:

   -  | **compute_new_temp**: if :math:`\rho e < 0`, then
      | call nyx_eos_given_RT (e,...,small_temp,...) in order to compute a new energy, :math:`e`.
      | This energy is then used to define a new :math:`E = e + ke`
      | Coming out of this the temperature is equal to small_temp and the energy :math:`e` has been reset.

   -  | **reset_internal_energy**: if :math:`e < 0` and :math:`E - ke < 0` then
      | call nyx_eos_given_RT (e,...,small_temp,...) in order to compute a new energy, :math:`e`. This energy is also used to define a new :math:`E = e + ke`

-  **large_temp**:

   -  | **compute_new_temp**: if :math:`T > \mathrm{large\_temp}`, and the input flag ``nyx.local_max_temp_dt=1`` then
      | set :math:`T = \mathrm{large\_temp}` and call nyx_eos_given_RT (e,...,large_temp,...) in order to compute a new energy, :math:`e`.
      | This energy is then used to define a new :math:`E = e + ke`
      | Coming out of this the temperature is equal to large_temp and the energy :math:`e` has been reset.

-  **small_pres**:

   -  | **tracexy / tracez / tracexy_ppm / tracez_ppm**:
      | qpres = max(qpres,small_pres)
      | for qxp/qyp, qyp/qym and qzp/qzm. This modifies the primitive density and pressure inside the tracing, not the underlying state variables

   -  **riemannus** – we set

      | pstar = max(small_pres,pstar)
      | pgdnv = max(small_pres,pgdnv). Note that pgdnv is the pressure explicitly used in the fluxes.

   -  **uflatten** – small_pres is used to keep the denominator away from zero

   -  Everywhere we define values of pressure on a face, we set that value to be at least small_pres.
