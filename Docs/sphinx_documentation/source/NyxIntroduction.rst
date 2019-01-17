**********************
Introduction to Nyx
**********************

Nyx is a adaptive mesh, hydrodynamics code that is
designed to model astrophysical reacting flows on massively parallel
computers.

The major capabilities:

  * 3-dimensional unsplit, 2nd-order hydrodynamics

  * adaptive mesh refinement with subcycling; jumps of 2x and 4x between levels

  * full Poisson gravity (with isolated boundary conditions)

  * parallelization via MPI + OpenMP

Units and Conventions
=====================

Nyx supports both CGS and Cosmological units. In the equation of state calls,
conversions must be made between CGS units and code units.
Table \ `[table:units] <#table:units>`__ shows some of the common symbols / names used
throughout the code documentation and papers.

.. table:: [table:units] Common quantities and units.

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

words

Inputs and Outputs
===================

We support two different systems of units in : CGS and Cosmological.
All inputs and problem initialization should be specified consistently with one of these sets of units.
No internal conversions of units occur within the code, so the output must be interpreted appropriately.
The default is cosmological units.
If you want to use CGS units instead, then set
USE_CGS = TRUE
in your GNUmakefile. This will select the file
constants_cgs.f90 instead of constants_cosmo.f90 from the
Nyx/constants directory.

.. table:: [table:inputs]
	   Units expected for inputs and outputs.
	   
   +-------------+-------------+-------------+-------------+-------------+
   |             |             |             |             |             |
   +-------------+-------------+-------------+-------------+-------------+
   | Location    | Variable    | CGS         | Cosmologica | Conversion  |
   |             |             |             | l           | Data        |
   +-------------+-------------+-------------+-------------+-------------+
   |             |             |             |             |             |
   +-------------+-------------+-------------+-------------+-------------+
   |             |             |             |             |             |
   +-------------+-------------+-------------+-------------+-------------+
   | inputs file | **geometry. | cm          | Mpc         | 1Mpc =      |
   |             | prob_lo**   |             |             | 3.08568025e |
   |             |             |             |             | 24          |
   |             |             |             |             | cm          |
   +-------------+-------------+-------------+-------------+-------------+
   |             | **geometry. | cm          | Mpc         | 1Mpc =      |
   |             | prob_hi**   |             |             | 3.08568025e |
   |             |             |             |             | 24          |
   |             |             |             |             | cm          |
   +-------------+-------------+-------------+-------------+-------------+
   |             |             |             |             |             |
   +-------------+-------------+-------------+-------------+-------------+
   |             |             |             |             |             |
   +-------------+-------------+-------------+-------------+-------------+
   | Hydro       | density     | g /         | M\ :math:`_ | 1           |
   | Initializat |             | cm\ :math:` | \odot`      | (M:math:`_\ |
   | ion         |             | ^3`         | /           | odot`       |
   |             |             |             | Mpc\ :math: | /           |
   |             |             |             | `^3`        | Mpc\ :math: |
   |             |             |             |             | `^3`)       |
   |             |             |             |             | =           |
   |             |             |             |             | .06769624e- |
   |             |             |             |             | 39          |
   |             |             |             |             | (g/cm:math: |
   |             |             |             |             | `^3`)       |
   +-------------+-------------+-------------+-------------+-------------+
   |             |             |             |             |             |
   +-------------+-------------+-------------+-------------+-------------+
   |             |             |             |             |             |
   +-------------+-------------+-------------+-------------+-------------+
   | Hydro       | velocities  | cm/s        | km/s        | 1km = 1.e5  |
   | Initializat |             |             |             | cm          |
   | ion         |             |             |             |             |
   +-------------+-------------+-------------+-------------+-------------+
   |             |             |             |             |             |
   +-------------+-------------+-------------+-------------+-------------+
   |             |             |             |             |             |
   +-------------+-------------+-------------+-------------+-------------+
   | Hydro       | momenta     | (g/cm:math: | (M:math:`_\ | 1km = 1.e5  |
   | Initializat |             | `^3`)       | odot`/Mpc:m | cm          |
   | ion         |             | (cm/s)      | ath:`^3`)   |             |
   |             |             |             | (km/s)      |             |
   +-------------+-------------+-------------+-------------+-------------+
   |             |             |             |             | 1           |
   |             |             |             |             | (M:math:`_\ |
   |             |             |             |             | odot`       |
   |             |             |             |             | /           |
   |             |             |             |             | Mpc\ :math: |
   |             |             |             |             | `^3`)       |
   |             |             |             |             | =           |
   |             |             |             |             | .06769624e- |
   |             |             |             |             | 39          |
   |             |             |             |             | g/cm\ :math |
   |             |             |             |             | :`^3`       |
   +-------------+-------------+-------------+-------------+-------------+
   |             |             |             |             |             |
   +-------------+-------------+-------------+-------------+-------------+
   |             |             |             |             |             |
   +-------------+-------------+-------------+-------------+-------------+
   | Hydro       | temperature | K           | K           | 1           |
   | Initializat |             |             |             |             |
   | ion         |             |             |             |             |
   +-------------+-------------+-------------+-------------+-------------+
   |             |             |             |             |             |
   +-------------+-------------+-------------+-------------+-------------+
   |             |             |             |             |             |
   +-------------+-------------+-------------+-------------+-------------+
   | Hydro       | specific    | erg/g=      | (km/s):math | 1           |
   | Initializat | energy      | (cm/s):math | :`^2`       | (km/s):math |
   | ion         | (:math:`e`  | :`^2`       |             | :`^2`       |
   |             | or          |             |             | = 1.e10     |
   |             | :math:`E`)  |             |             | (cm/s):math |
   |             |             |             |             | :`^2`       |
   +-------------+-------------+-------------+-------------+-------------+
   |             |             |             |             |             |
   +-------------+-------------+-------------+-------------+-------------+
   |             |             |             |             |             |
   +-------------+-------------+-------------+-------------+-------------+
   | Hydro       | energy      | erg /       | (M:math:`_\ | 1           |
   | Initializat | (:math:`\rh | cm\ :math:` | odot`/Mpc:m | (km/s):math |
   | ion         | o e`        | ^3 =`       | ath:`^3`)   | :`^2`       |
   |             | or          |             | (km/s):math | = 1.e10     |
   |             | :math:`\rho |             | :`^2`       | (cm/s):math |
   |             |  E`)        |             |             | :`^2`       |
   +-------------+-------------+-------------+-------------+-------------+
   |             |             | (g/cm:math: |             | 1           |
   |             |             | `^3`)       |             | (M:math:`_\ |
   |             |             | (cm/s):math |             | odot`       |
   |             |             | :`^2`       |             | /           |
   |             |             |             |             | Mpc\ :math: |
   |             |             |             |             | `^3`)       |
   |             |             |             |             | =           |
   |             |             |             |             | .06769624e- |
   |             |             |             |             | 39          |
   |             |             |             |             | g/cm\ :math |
   |             |             |             |             | :`^3`       |
   +-------------+-------------+-------------+-------------+-------------+
   |             |             |             |             |             |
   +-------------+-------------+-------------+-------------+-------------+
   |             |             |             |             |             |
   +-------------+-------------+-------------+-------------+-------------+
   | Particle    | particle    | g           | M\ :math:`_ | 1           |
   | Initializat | mass        |             | \odot`      | M\ :math:`_ |
   | ion         |             |             |             | \odot`      |
   |             |             |             |             | =           |
   |             |             |             |             | 1.98892e33  |
   |             |             |             |             | g           |
   +-------------+-------------+-------------+-------------+-------------+
   |             |             |             |             |             |
   +-------------+-------------+-------------+-------------+-------------+
   |             |             |             |             |             |
   +-------------+-------------+-------------+-------------+-------------+
   | Particle    | particle    | cm          | Mpc         | 1 Mpc =     |
   | Initializat | locations   |             |             | 3.08568025e |
   | ion         |             |             |             | 24          |
   |             |             |             |             | cm          |
   +-------------+-------------+-------------+-------------+-------------+
   |             |             |             |             |             |
   +-------------+-------------+-------------+-------------+-------------+
   |             |             |             |             |             |
   +-------------+-------------+-------------+-------------+-------------+
   | Particle    | particle    | cm/s        | km/s        | 1 km = 1e5  |
   | Initializat | velocities  |             |             | cm          |
   | ion         |             |             |             |             |
   +-------------+-------------+-------------+-------------+-------------+
   |             |             |             |             |             |
   +-------------+-------------+-------------+-------------+-------------+
   |             |             |             |             |             |
   +-------------+-------------+-------------+-------------+-------------+
   | Output      | Pressure    | g           | M\ :math:`_ | 1           |
   |             |             | (cm/s):math | \odot`      | M\ :math:`_ |
   |             |             | :`^2`       | (km/s):math | \odot`      |
   |             |             | /           | :`^2`       | (km/s):math |
   |             |             | cm\ :math:` | /           | :`^2`       |
   |             |             | ^3`         | Mpc\ :math: | /           |
   |             |             |             | `^3`        | Mpc\ :math: |
   |             |             |             |             | `^3`        |
   |             |             |             |             | =           |
   +-------------+-------------+-------------+-------------+-------------+
   |             |             |             |             | .06769624e- |
   |             |             |             |             | 29          |
   |             |             |             |             | g           |
   |             |             |             |             | (cm/s):math |
   |             |             |             |             | :`^2`       |
   |             |             |             |             | /           |
   |             |             |             |             | cm\ :math:` |
   |             |             |             |             | ^3`         |
   +-------------+-------------+-------------+-------------+-------------+
   |             |             |             |             |             |
   +-------------+-------------+-------------+-------------+-------------+
   |             |             |             |             |             |
   +-------------+-------------+-------------+-------------+-------------+
   | Output      | Gravity     | (cm/s) / s  | (km/s):math | 1           |
   |             |             |             | :`^2`       | M\ :math:`_ |
   |             |             |             | / Mpc       | \odot`      |
   |             |             |             |             | (km/s):math |
   |             |             |             |             | :`^2`       |
   |             |             |             |             | /           |
   |             |             |             |             | Mpc\ :math: |
   |             |             |             |             | `^3`        |
   |             |             |             |             | =           |
   +-------------+-------------+-------------+-------------+-------------+
   |             |             |             |             |             |
   +-------------+-------------+-------------+-------------+-------------+
   |             |             |             |             |             |
   +-------------+-------------+-------------+-------------+-------------+
   | Output      | Time        | s           | (Mpc/km) s  | 1 Mpc =     |
   |             |             |             |             | 3.08568025e |
   |             |             |             |             | 19          |
   |             |             |             |             | km          |
   +-------------+-------------+-------------+-------------+-------------+
   |             |             |             |             |             |
   +-------------+-------------+-------------+-------------+-------------+

[Table:Inputs]

.. table:: [table:constants]
	   Physical constant values
   
   +-----------------+-----------------+-----------------+-----------------+
   |                 |                 |                 |                 |
   +-----------------+-----------------+-----------------+-----------------+
   | Constant        | CGS             | Cosmological    | Conversion Data |
   +-----------------+-----------------+-----------------+-----------------+
   |                 |                 |                 |                 |
   +-----------------+-----------------+-----------------+-----------------+
   |                 |                 |                 |                 |
   +-----------------+-----------------+-----------------+-----------------+
   | Gravitational   | 6.67428e-8 cm   | 4.3019425e-9        |                 |
   | constant        | (cm/s)          | Mpc                 |                 |
   |                 | :math:`^2` g\   | (km/s)              |                 |
   | (:math:`G`)     |                 | :math:`^2`          |                 |
   |                 | :math:`^{-1}`   | M\                  |                 |
   |                 |                 | :math:`_\odot^{-1}` |                 |
   |                 |                 |                     |                 |
   +-----------------+-----------------+-----------------+-----------------+
   |                 |                 |                 |                 |
   +-----------------+-----------------+-----------------+-----------------+
   |                 |                 |                 |                 |
   +-----------------+-----------------+-----------------+-----------------+
   | Avogadro’s      | 6.02214129e23   | 1.1977558e57    | 1               |
   | number          | g\ :math:`^{-1}`| M\ :math:`_\odo | M\ :math:`_\odo |
   | (:math:`n_A`)   |                 | t^{-1}`         | t`              |
   |                 |                 |                 | = 1.98892e33 g  |
   +-----------------+-----------------+-----------------+-----------------+
   |                 |                 |                 |                 |
   +-----------------+-----------------+-----------------+-----------------+
   |                 |                 |                 |                 |
   +-----------------+-----------------+-----------------+-----------------+
   | Boltzmann’s     | 1.3806488e-16   | 0.6941701e-59   | 1               |
   | constant        | erg / K         | M\ :math:`_\odo | M\ :math:`_\odo |
   | (:math:`k_B`)   |                 | t`              | t`              |
   |                 |                 | (km/s):math:`^2 | (km/s):math:`^2 |
   |                 |                 | `               | `               |
   |                 |                 | / K             | = 1.98892e43 g  |
   |                 |                 |                 | (cm/s):math:`^2 |
   |                 |                 |                 | `               |
   +-----------------+-----------------+-----------------+-----------------+
   |                 |                 |                 |                 |
   +-----------------+-----------------+-----------------+-----------------+
   |                 |                 |                 |                 |
   +-----------------+-----------------+-----------------+-----------------+
   | Hubble constant | 100 (km/s) /    | 32.407764868e-1 | 1 Mpc =         |
   | (:math:`H`)     | Mpc             | 9               | 3.08568025e19   |
   |                 |                 | s\ :math:`^{-1} | km              |
   |                 |                 | `               |                 |
   +-----------------+-----------------+-----------------+-----------------+
   |                 |                 |                 |                 |
   +-----------------+-----------------+-----------------+-----------------+

.. raw:: latex

   \clearpage

The only other place that dimensional numbers are used in the code is in the tracing and Riemann solve.
We set three *small* numbers which need to be consistent with the data specified.
Each of these can be specified in the inputs file.

-  small_dens – small value for density

-  small_p – small value for pressure

-  small_T – small value for temperature

These are the places that each is used in the code:

-  **small_dens**

   -  | **subroutine enforce_minimum_density** (called after subroutine consup) – if :math:`\rho <` small_dens then :math:`\rho` is set to the
        minimum value of the 26 neighbors. This also modifies momenta, :math:`\rho E` and :math:`\rho e` so that velocties, :math:`E` and :math:`e` remain unchanged.

   -  | **subroutine tracexy / tracez / tracexy_ppm / tracez_ppm**:
      | qxp = max(qxp,small_dens)
      | qxm = max(qxm,small_dens)
      | and analogously for qyp/qym and qzp/qzm. This only modifies density inside the tracing, not the other variables

   -  **subroutine riemannus** – we set

      wsmall = small_dens \* csmall

      and then

      | wl = max(wsmall, sqrt(gaml \* pl \* rl))
      | wr = max(wsmall, sqrt(gamr \* pr \* rr))

      Also, we set

      ro = max(small_dens,ro)

      where ro = 0.5 \* (rl + rr) – this state is only chosen when ustar = 0, and

      rstar = max(small_dens,rstar)

      where rstar = ro + (pstar-po)/co:math:`^2`

   -  **subroutine react_state** – only compute reaction if :math:`\rho >` small_dens

-  **small_temp**:

   -  | **subroutine ctoprim**: if :math:`\rho e < 0`, then
      | call subroutine nyx_eos_given_RTX (e,...,small_temp,...) in order to compute a new energy, :math:`e`.
      | This energy is then used to
      | call subroutine nyx_eos_given_ReX in order to compute the sound speed, :math:`c.`
      | Coming out of this the temperature is equal to small_temp and the energy :math:`e` has been reset.

   -  | **subroutine react_state**: if :math:`\rho e < 0`, then
      | call subroutine nyx_eos_given_RTX (e,...,small_temp,...) in order to compute a new energy, :math:`e`.
      | This energy is then used to proceed with the burner routine.

   -  | **subroutine reset_internal_energy**: if :math:`e < 0` and :math:`E - ke < 0` then
      | call subroutine nyx_eos_given_RTX (e,...,small_temp,...) in order to compute a new energy, :math:`e`. This energy is also used to
      | define a new :math:`E = e + ke`

-  **small_pres**:

   -  **subroutine riemannus** – we set

      | pstar = max(small_pres,pstar)
      | pgdnv = max(small_pres,pgdnv). Note that pgdnv is the pressure explicitly used in the fluxes.

   -  **subroutine uflaten** – small_pres is used to keep the denominator away from zero

   -  Everywhere we define values of pressure on a face, we set that value to be at least small_pres.
