******************
Stochastic Forcing
******************

There is an option to apply a stochastic force field. 

See Nyx/Exec/DrivenTurbulence for an example; note that ::

    nyx.do_forcing = 1

must be set in the inputs file.

The external forcing term in the momentum equation 
(`[eq:momt] <#eq:momt>`__) is then given by

  .. math:: {\bf S}_{\rho \Ub} = \rho_b \mathbf{f}

  where the acceleration field :math:`\mathbf{f}(\mathbf{x},t)` is
  computed as inverse Fourier transform of the forcing spectrum
  :math:`\widehat{\mathbf{f}}(\mathbf{k},t`). The time evolution of each
  wave mode is given by an Ornstein-Uhlenbeck process (see
  :raw-latex:`\cite{SchmHille06,Schmidt14}` for details). Since the real
  space forcing acts on large scales :math:`L`, non-zero modes are
  confined to a narrow window of small wave numbers with a prescribed
  shape (the forcing profile). The resulting flow reaches a
  statistically stationary and isotropic state with a root-mean-square
  velocity of the order :math:`V=L/T`, where the integral time scale
  :math:`T` (also known as large-eddy turn-over time) is usually set
  equal to the autocorrelation time of the forcing. It is possible to
  vary the force field from solenoidal (divergence-free) if the weight
  parameter :math:`\zeta=1` to dilational (rotation-free) if
  :math:`\zeta=0`.

To maintain a nearly constant root-mean-square Mach number, a simple
model for radiative heating and cooling around a given equilibrium
temperature :math:`T_0` is applied in the energy
equation (`[eq:energy] <#eq:energy>`__):

.. math:: S_{\rho E} = S_{\rho e} + \Ub \cdot {\bf S}_{\rho \Ub} = -\frac{\alpha k_{\rm B}(T-T_0)}{\mu m_{\rm H}(\gamma-1)} + \rho_b\Ub\cdot\mathbf{f}

The parameters :math:`T_0` and :math:`\alpha` correspond to temp0 and
alpha, respectively, in the probin file (along with rho0 for the mean
density, which is unity by default). While the gas is adiabatic for
:math:`\alpha=0`, it becomes nearly isothermal if the cooling time scale
given by :math:`1/\alpha` is chosen sufficiently short compared to
:math:`T`. For performance reasons, a constant composition
(corresponding to constant molecular weight :math:`\mu`) is assumed.

List of Parameters
==================

+-------------------------+--------------------+-----------------+-------------+
| Parameter               | Definition         | Acceptable      | Default     |
|                         |                    | Values          |             |
+=========================+====================+=================+=============+
| **forcing.seed**        | seed of the        | Integer         | 27011974    |
|                         | random number      | :math:`>0`      |             |
|                         | generator          |                 |             |
+-------------------------+--------------------+-----------------+-------------+
| **forcing.profile**     | shape of           | 1 (plane), 2    | 3           |
|                         | forcing            | (band), 3       |             |
|                         | spectrum           | (parabolic)     |             |
+-------------------------+--------------------+-----------------+-------------+
| **forcing.alpha**       | ratio of domain    | Integer         | 2 2 2       |
|                         | size :math:`X`     | :math:`>0`      |             |
|                         | to integral        |                 |             |
|                         | length             |                 |             |
|                         | :math:`L=X/\alpha` |                 |             |
|                         |                    |                 |             |
+-------------------------+--------------------+-----------------+-------------+
| **forcing.band_width**  | band width of      | Real            | 1.0 1.0 1.0 |
|                         | the forcing        | :math:`\ge 0`   |             |
|                         | spectrum           | and             |             |
|                         | relative to        | :math:`\le 1`   |             |
|                         | alpha              |                 |             |
+-------------------------+--------------------+-----------------+-------------+
| **forcing.intgr_vel**   | characteristic     | Real            | must be set |
|                         | velocity           | :math:`> 0`     |             |
|                         | :math:`V`          |                 |             |
+-------------------------+--------------------+-----------------+-------------+
| **forcing.auto_corrl**  | autocorrelation    | Real            | 1.0 1.0 1.0 |
|                         | time in units      | :math:`> 0`     |             |
|                         | of                 |                 |             |
|                         | :math:`T=L/V`      |                 |             |
+-------------------------+--------------------+-----------------+-------------+
| **forcing.soln_weight** | weight             | Real            | 1.0         |
|                         | :math:`\zeta`      | :math:`\ge 0`   |             |
|                         | of solenoidal      | and             |             |
|                         | relative to        | :math:`\le 1`   |             |
|                         | dilatational       |                 |             |
|                         | modes              |                 |             |
+-------------------------+--------------------+-----------------+-------------+

Triples for forcing.alpha, forcing.band_width, forcing.intgr_vel, and
forcing.auto_corrl correspond to the three spatial dimensions.
