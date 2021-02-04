
Units and Conventions
=====================

Nyx supports both CGS and Cosmological units. 

All inputs and problem initialization should be specified consistently with one of these sets of units.
No internal conversions of units occur within the code, so the output must be interpreted appropriately.
The default is cosmological units.

If you want to use CGS units instead, then set

**USE_CGS = TRUE**

in your GNUmakefile. This will select the file constants_cgs.f90 instead of constants_cosmo.f90 from the
Nyx/constants directory.

In the equation of state calls,
conversions must be made between CGS units and code units.
:numref:`table:units` shows some of the common symbols / names used
throughout the code documentation and papers.

.. _table:units:
.. table:: Common quantities and units.

   +-----------------------+-----------------------+-----------------------+
   | name                  | units                 | description           |
   +=======================+=======================+=======================+
   | :math:`t`             | s                     | time                  |
   +-----------------------+-----------------------+-----------------------+
   | :math:`\rho`          | :math:`\gcc`          | mass density          |
   +-----------------------+-----------------------+-----------------------+
   | :math:`\ub`           | :math:`\cms`          | velocity vector       |
   +-----------------------+-----------------------+-----------------------+
   | :math:`p`             | :math:`\presunit`     | pressure              |
   +-----------------------+-----------------------+-----------------------+
   | :math:`\gb`           | :math:`\accelunit`    | gravitational         |
   |                       |                       | acceleration          |
   +-----------------------+-----------------------+-----------------------+
   | :math:`\Sb`           | varies                | source term           |
   +-----------------------+-----------------------+-----------------------+
   | :math:`S_{\Lambda}`   | varies                | source term           |
   +-----------------------+-----------------------+-----------------------+
   | :math:`E`             | :math:`\ergg`         | specific total energy |
   +-----------------------+-----------------------+-----------------------+
   | :math:`e`             | :math:`\ergg`         | specific internal     |
   |                       |                       | energy                |
   +-----------------------+-----------------------+-----------------------+
   | :math:`T`             | :math:`K`             | temperature           |
   +-----------------------+-----------------------+-----------------------+

:numref:`table:inputs` associates the numerical input data with the corresponding unit scaling.

.. _table:inputs:
.. table:: 
	   Units expected for inputs and outputs.
	   
   ======================= ========================================= ===================================== ====================================================== ==========================================================================
   Location                Variable                                  CGS                                   Cosmological                                           Conversion Data
   ======================= ========================================= ===================================== ====================================================== ==========================================================================
   inputs file             **geometry.prob_lo**                      cm                                    Mpc                                                    1Mpc = 3.08568025e24 cm
   \                       **geometry.prob_hi**                      cm                                    Mpc                                                    1Mpc = 3.08568025e24 cm
   Hydro Initialization    density                                   g / cm\ :math:`^3`                    M\ :math:`_\odot` / Mpc\ :math:`^3`                    1 (M\ :math:`_\odot` / Mpc\ :math:`^3`) = .06769624e-39 (g/cm\ :math:`^3`)
   Hydro Initialization    velocities                                cm/s                                  km/s                                                   1km = 1.e5 cm
   Hydro Initialization    momenta                                   (g/cm\ :math:`^3`) (cm/s)             (M\ :math:`_\odot`/Mpc\ :math:`^3`) (km/s)             1km = 1.e5 cm
   \                                                                                                                                                              1 (M\ :math:`_\odot` / Mpc\ :math:`^3`) = .06769624e-39 g/cm\ :math:`^3`
   Hydro Initialization    temperature                               K                                     K                                                      1
   Hydro Initialization    specific energy (:math:`e` or :math:`E`)  erg/g= (cm/s)\ :math:`^2`             (km/s)\ :math:`^2`                                     1 (km/s)\ :math:`^2` = 1.e10 (cm/s)\ :math:`^2`
   Hydro Initialization    energy (:math:`\rho e` or :math:`\rho E`) erg / cm\ :math:`^3 =`                (M\ :math:`_\odot`/Mpc\ :math:`^3`) (km/s)\ :math:`^2` 1 (km/s)\ :math:`^2` = 1.e10 (cm/s)\ :math:`^2`
   \                                                                 (g/cm\ :math:`^3`) (cm/s)\ :math:`^2`                                                        1 (M\ :math:`_\odot` / Mpc\ :math:`^3`) = .06769624e-39 g/cm\ :math:`^3`
   Particle Initialization particle mass                             g                                     M\ :math:`_\odot`                                      1 M\ :math:`_\odot` = 1.98892e33 g
   Particle Initialization particle locations                        cm                                    Mpc                                                    1 Mpc = 3.08568025e24 cm
   Particle Initialization particle velocities                       cm/s                                  km/s                                                   1 km = 1e5 cm
   Output                  Pressure                                  g (cm/s)\ :math:`^2` / cm\ :math:`^3` M\ :math:`_\odot` (km/s)\ :math:`^2` / Mpc\ :math:`^3` 1 M\ :math:`_\odot` (km/s)\ :math:`^2` / Mpc\ :math:`^3` =
   \                                                                                                                                                              .06769624e-29 g (cm/s)\ :math:`^2` / cm\ :math:`^3`
   Output                  Gravity                                   (cm/s) / s                            (km/s)\ :math:`^2` / Mpc                               1 M\ :math:`_\odot` (km/s)\ :math:`^2` / Mpc\ :math:`^3` =
   Output                  Time                                      s                                     (Mpc/km) s                                             1 Mpc = 3.08568025e19 km
   ======================= ========================================= ===================================== ====================================================== ==========================================================================

:numref:`table:constants` lists the values used for cosmological constants in the code with their associated units. Note that :math:`\Omega_m`, :math:`\Omega_b`, :math:`\Omega_r`  and :math:`h` are set in the inputs file.
   
.. _table:constants:
.. table::
	   Physical constant values

   ================================== ================================================= ========================================================== ========================================================================
   Constant                           CGS                                               Cosmological                                               Conversion Data
   ================================== ================================================= ========================================================== ========================================================================
   Gravitational constant (:math:`G`) 6.67428e-8 cm (cm/s)\ :math:`^2` g\ :math:`^{-1}` 4.3019425e-9 Mpc (km/s)\ :math:`^2` M\ :math:`_\odot^{-1}` 
   Avogadro’s number (:math:`n_A`)    6.02214129e23 g\ :math:`^{-1}`                    1.1977558e57 M\ :math:`_\odot^{-1}`                        1 M\ :math:`_\odot` = 1.98892e33 g
   Boltzmann’s constant (:math:`k_B`) 1.3806488e-16 erg / K                             0.6941701e-59 M\ :math:`_\odot` (km/s)\ :math:`^2` / K     1 M\ :math:`_\odot` (km/s)\ :math:`^2` = 1.98892e43 g (cm/s)\ :math:`^2`
   Hubble constant (:math:`H`)        100 (km/s) / Mpc                                  32.407764868e-19 s\ :math:`^{-1}`                          1 Mpc = 3.08568025e19 km
   ================================== ================================================= ========================================================== ========================================================================

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
