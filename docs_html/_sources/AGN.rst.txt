The AGN Model in Nyx
====================

In the AGN model, super-massive black hole (SMBH) particles are formed
at *haloes*, where each halo is defined by a connected mass enclosed by
a user-defined density isocontour. In order to find haloes, we use the
Reeber package described in Section \ `2 <#sec:Reeber>`__. Each AGN
particle has the standard dark matter particle attributes of position,
velocity, and mass, as well as two additional attributes, its stored
accretion energy and its mass accretion rate.

.. table:: Parameters of the AGN model

   ================ ============================ ======================== ===========================================
   In “probin” file Parameter                    Fiducial value           Explanation
   ================ ============================ ======================== ===========================================
   \*               :math:`M_{\rm h, min}`       :math:`10^{10}~ M_\odot` Minimum halo mass for SMBH placement
   \*               :math:`M_{\rm seed}`         :math:`10^5~ M_\odot`    Seed mass of SMBH
   T_min            :math:`T_{\rm min}`          :math:`10^7` K           Minimum heating of the surrounding gas
   bondi_boost      :math:`\alpha`               100                      Bondi accretion boost factor
   max_frac_removed :math:`f_{\rm max, removed}` 0.5                      Maximum fraction of mass removed from gas
   eps_rad          :math:`\epsilon_{\rm r}`     0.1                      Radiation efficiency
   eps_coupling     :math:`\epsilon{\rm c}`      0.15                     Coupling efficiency
   eps_kinetic      :math:`\epsilon_{\rm kin}`   0.1                      Kinetic feedback efficiency
   frac_kinetic     :math:`f_{\rm kin}`          0                        Fraction of feedback energy that is kinetic
   ================ ============================ ======================== ===========================================

[tab:agn_params1]

\* :math:`M_{\rm h, min}` and :math:`M_{\rm seed}` are not set in the
“probin” file, but in the inputs file, by respectively Nyx.mass_halo_min
and Nyx.mass_seed.

Creating AGN Particles from Haloes
----------------------------------

Each halo with threshold mass of :math:`M_h \geqslant M_{\rm h, min}`
that does not already host a black hole particle is seeded with a black
hole of mass :math:`M_{\rm seed}`. The initial position of this AGN
particle is the center of the cell where the density is highest in the
halo.

When an AGN particle is created, the density in its cell is reduced by
the amount required for mass to be conserved, and the velocity of the
AGN particle is initialized so that momentum is conserved. The accretion
energy and mass accretion rate are initialized to zero.

Merging AGN Particles
---------------------

Two AGN particles merge when both of these conditions obtain:

#. The distance between them, :math:`l`, is less than the mesh spacing,
   :math:`h`.

#. [velocity-item] The difference of their velocities,
   :math:`v_{\rm rel}`, is less than the circular velocity at distance
   :math:`l`:

   .. math:: v_{\rm rel} < \sqrt{GM_{\rm BH}/l}

   where :math:`M_{\rm BH}` is the mass of the more massive SMBH in the
   pair, and :math:`G` is the gravitational constant.

Criterion \ `[velocity-item] <#velocity-item>`__ above is necessary in
order to prevent AGN particles from merging during a fly-through
encounter of two haloes, as this could lead to AGN particles being
quickly removed from the host halo due to momentum conservation.

The merger of two AGN particles is implemented as the less massive one
being removed, and its mass and momentum being transferred to the more
massive one.

Accretion
---------

For an AGN particle of mass :math:`M_{\rm BH}`, the Bondi–Hoyle
accretion rate is

.. math::

   \dot{M}_{\rm B} = \alpha
   \frac{4 \pi G^2 M_{\rm BH}^2 \overline{\rho}}{(\overline{c_s^2} + \overline{u^2})^{3/2}} ,

where :math:`\overline{\rho}`, :math:`\overline{c_s^2}`, and
:math:`\overline{u^2}` are volume averages with a cloud-in-cell stencil
of the gas’s density, squared sound speed, and squared velocity,
respectively, in the neighborhood of the particle.

The maximum black hole accretion rate is the Eddington limit,

.. math::

   \dot{M}_{\rm Edd} = 
   \frac{4 \pi G M_{\rm BH} m_{\rm p}}{\epsilon_{\rm r} \sigma_{\rm T} c} \, ,

with proton mass :math:`m_{\rm p}`, Thomson cross section
:math:`\sigma_{\rm T}`, and speed of light :math:`c`.

The mass accretion rate of the SMBH is the smaller of the two rates
above:
:math:`\dot{M}_{\rm acc} = {\rm min} \{ \dot{M}_{\rm B}, \dot{M}_{\rm Edd} \}`.
Then the gas will lose mass :math:`\dot{M}_{\rm acc} \Delta t`, where
:math:`\Delta t` is the length of the time step. However,
:math:`\dot{M}_{\rm acc}` is adjusted downward if necessary so that when
cloud-in-cell stencil weights are applied in the neighborhood of the
particle, the fraction of gas density removed from any cell of the
stencil is at most :math:`f_{\rm max, removed}`.

The mass of the AGN particle increases by
:math:`(1-\epsilon_{\rm r}) \dot{M}_{\rm acc} \Delta t`, while
:math:`\dot{M}_{\rm acc} \Delta t` amount of gas mass is removed from
the grid according to cloud-in-cell stencil weights in the neighborhood
of the particle. The momentum transfer can be computed by assuming the
velocity of the gas is unchanged; thus the gas in each cell loses
momentum in proportion to its mass loss, and the particle gains the sum
of the gas momentum loss multiplied by :math:`(1-\epsilon_{\rm r})`.

Feedback Energy
---------------

Feedback energy is stored in an AGN particle variable
:math:`E_{\rm AGN}`, and is accumulated over time until released. The
fraction :math:`f_{\rm kin}` goes to kinetic energy, and the rest to
thermal energy.

Thermal Feedback
~~~~~~~~~~~~~~~~

We increment :math:`E_{\rm AGN}` by thermal feedback energy, calculated
from the mass accretion rate as

.. math::

   % high energy
   \Delta E_{\rm thermal} = (1 - f_{\rm kin})
   \epsilon_{\rm c} \epsilon_{\rm r} \dot{M}_{\rm acc} c^2 \Delta t .

Kinetic/Momentum Feedback
~~~~~~~~~~~~~~~~~~~~~~~~~

We increment :math:`E_{\rm AGN}` by the kinetic feedback energy

.. math::

   % low energy
   \Delta E_{\rm kinetic} =
   f_{\rm kin} \epsilon_{\rm kin} \dot{M}_{\rm acc} c^2 \Delta t .

We also need to adjust the energy density and momentum density of the
gas. We do this by computing a jet velocity

.. math:: \vec{v}_{\rm jet} = \sqrt{\frac{2 \Delta E_{\rm kinetic}}{m_{\rm g}}} \vec{n}

where :math:`m_{\rm g}` is the total gas mass inside the cloud-in-cell
local environment, and :math:`\vec{n}` is a *randomly* chosen unit
vector. We add :math:`\rho \vec{v}` to the momentum density
:math:`\vec{p}` of the gas, and :math:`\vec{v}_{\rm jet} \cdot \vec{p}`
to its energy density, both of these weighted by the cloud-in-cell
stencil of the particle.

Releasing Feedback Energy
~~~~~~~~~~~~~~~~~~~~~~~~~

The accumulated energy is released when

.. math::

   E_{\rm AGN} > m_{\rm g} \overline{e}
   \label{eq:E_agn}

where :math:`\overline{e}` is the average specific internal energy of
the gas over the cloud-in-cell stencil, obtained from the equation of
state using temperature :math:`T_{\rm min}` and average density of the
gas over the same stencil, and :math:`m_{\rm g}` is the total gas mass
inside the cloud-in-cell local environment.

.. _sec:Reeber:

The Reeber Package
==================

Reeber is a separate package with a halo finder. Here are the Reeber
parameters that are assigned in the input file.

=================================== =================================== ===================================== =========
Parameter                           Definition                          Acceptable Values                     Default
=================================== =================================== ===================================== =========
**reeber.halo_int**                 timesteps between halo finder calls Integer                               -1 (none)
**reeber.negate**                   allow negative values for analysis  0 if false, 1 if true                 1
**reeber.halo_density_vars**        density variable list               density, particle_mass_density        “density”
**reeber.halo_extrema_threshold**   extrema threshold for haloes        Real                                  200.
**reeber.halo_component_threshold** component threshold for haloes      Real                                  82.
**reeber.absolute_halo_thresholds** are halo thresholds absolute        0 if multiples of mean, 1 if absolute 0
=================================== =================================== ===================================== =========

[Table:Reeber-inputs]
