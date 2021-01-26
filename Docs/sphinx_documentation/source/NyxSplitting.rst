============================================
Hydrodynamical and Heating-Cooling Splitting
============================================

We solve the equations of gas dynamics in a coordinate system that is comoving
with the expanding universe, with expansion factor, :math:`a,` related to the redshift, :math:`z`, by :math:`a = 1 / (1 + z).`

We describe the state of the gas
as :math:`\overline{{\bf U}} = (\rho, a \rho {\bf U}, a^2 \rho E, a^2 \rho e),`
then write the evolution of the gas as

.. math:: \frac{\partial\overline{{\bf U}}}{\partial t} = -\nabla\cdot{\bf F}+ S_e + S_g + S_{HC},

where :math:`{\bf F}= (1/a \; \rho {\bf U}, \rho {\bf U}{\bf U}, a (\rho {\bf U}E + p {\bf U}), a \rho {\bf U}e)`
is the flux vector,
:math:`S_e = (0, 0, 0, -a p \nabla \cdot {\bf U})` represents the additional term in the evolution
equation for internal energy, :math:`S_g = (0, \rho_b {\bf g}, a \rho {\bf U}\cdot {\bf g}, 0)`
represents the gravitational source terms,
and :math:`S_{HC} = (0, 0, a \rho \Lambda_{HC}, a \rho \Lambda_{HC})`
represents the combined heating and cooling source terms. The state, :math:`\overline{{\bf U}},` and
all source terms are defined at cell centers; the fluxes are defined on cell faces.

We compute :math:`{\bf F}`
using an unsplit Godunov method with characteristic tracing and full
corner coupling.

We track derived variables such as temperature and electron fraction ne based on the internal energy. The concentrations
of different isotopes of hydrogen and helium then depend on a Newton solve of the heating/cooling equations of state, which
are closely tied to the internal energy.

Strang Splitting
----------------

The original splitting used in Nyx is Strang splitting, where a half-step of the heating-cooling
is evolved, then a full step of the hydrodynamical terms, followed by a half-step of the heating-cooling.
This algorithm is the classical Strang splitting, adapted to include gravity and other forcing terms.

.. math::

   \label{eq:dens}
   \frac{\partial \rho}{\partial t} = - \frac{1}{a} \nabla \cdot (\rho {\bf U}) , \\

.. math::

   \begin{aligned}
   \label{eq:momt}
   \frac{\partial (a \rho {\bf U})}{\partial t} &=& 
   -             \nabla \cdot (\rho {\bf U}{\bf U}) 
   -             \nabla p 
   +             \rho {\bf g}, \end{aligned}

.. math::

   \begin{aligned}
   \frac{\partial (a^2 \rho E)}{\partial t} &=&  a \left[
    -\nabla \cdot (\rho {\bf U}E + p {\bf U})
   +  \rho {\bf U}\cdot {\bf g}
   +  \rho \Lambda_{HC}  \right]  . \end{aligned}

.. math::

   \begin{aligned}
   \frac{\partial (a^2 \rho e)}{\partial t} &=& a \left[ 
   - \nabla \cdot (\rho {\bf U}e)
   -  p \nabla \cdot {\bf U}
   +  \rho \Lambda_{HC}  \right]  . \end{aligned}

The algorithm at a single level of refinement begins by computing the time step, :math:`\Delta t,`
and advancing :math:`a` from :math:`t^n` to :math:`t^{n+1} = t^n + \Delta t`. The rest of the time step is
then composed of the following steps:

Step 1:
   Compute :math:`{\phi}^n` and :math:`{\bf g}^n` using :math:`\rho_b^n` and :math:`\rho_{dm}^n`,
   where :math:`\rho_{dm}^{n}` is computed from the particles at :math:`{\bf x}_i^{n}`.

   We note that in the single-level algorithm we can instead use :math:`{\bf g}` as computed at the
   end of the previous step because there have been no changes to :math:`{\bf x}_i`
   or :math:`\rho` since then.

Step 2:
   Interpolate :math:`{\bf g}^n` from the grid to the particle locations, then
   advance the particle velocities by :math:`\Delta t/ 2` and particle positions by :math:`\Delta t`.

   .. math::

      \begin{aligned}
           {\bf u}_i^{{n+\frac{1}{2}}} &=& \frac{1}{a^{{n+\frac{1}{2}}}} ((a^n {\bf u}^n_i) + \frac{\Delta t}{2} \; {\bf g}^n_i) \\
           {\bf x}_i^{n+1}  &=& {\bf x}^n_i + \frac{\Delta t}{a^{{n+\frac{1}{2}}}} {\bf u}_i^{{n+\frac{1}{2}}}\end{aligned}

Step 3:
   Advance :math:`{\bf U}` by :math:`\frac{1}{2}\Delta t` for first strang

   We advance :math:`e` by integrating the source terms in time for :math:`\frac{1}{2}\Delta t`

   .. math::

      \begin{aligned}
           ( e)^{n,\ast} &=& ( e)^n +  \int \Lambda_{HC} \; dt^\prime  .\end{aligned}

   We update :math:`(\rho e)^{n,\ast}=(\rho e)^{n,\ast}+\rho^{n}\left((e)^{n,\ast}-(e)^{n}\right)`

   We update :math:`(\rho E)^{n,\ast}=(\rho E)^{n,\ast}+\rho^{n}\left((e)^{n,\ast}-(e)^{n}\right)`

Step 4:
   Advance :math:`{\bf U}` by :math:`\Delta t` for advective terms

   Advance the solution using time-centered fluxes and :math:`S_e`
   and an explicit representation of :math:`S_g` at time :math:`t^n`:

   .. math:: {\bf U}^{n+1,\ast} = {\bf U}^{n,\ast} + \Delta tA^{n+\frac{1}{2}}+ \Delta tS_g^n

   where :math:`A^{n+\frac{1}{2}}` is computed by predicting from the :math:`{\bf U}^{n,\ast}` states.

   After adding the advective update terms :math:`A^{n+\frac{1}{2}}` and before calculating :math:`S_{g}`, apply a correction to :math:`{\bf U}^{n+1,\ast}` to enforce :math:`\rho > 1.1 \times small_dense`
		 
Step 5: 
   Second Strang step

   We advance :math:`e` by integrating the source terms in time for :math:`\frac{1}{2}\Delta t`

   .. math::

        \begin{aligned}
        ( e)^{n+1} &=& ( e)^{n+1,\ast } +  \int \Lambda_{HC} \; dt^\prime .\end{aligned}

   We update :math:`(\rho e)^{n+1}=(\rho e)^{n+1,\ast}+\rho^{n+1}\left((e)^{n+1}-(e)^{n+1,\ast}\right)`

   We update :math:`(\rho E)^{n+1}=(\rho E)^{n+1,\ast}+\rho^{n+1}\left((e)^{n+1}-(e)^{n+1,\ast}\right)`

   We store Ne and Temp based on eos\_ hc updates from :math:`(e)^{n+1}`

Step 6:
   Compute :math:`{\phi}^{n+1}` and :math:`{\bf g}^{n+1}` using
   :math:`\rho^{n+1,*}` and :math:`\rho_{dm}^{n+1}`, where :math:`\rho_{dm}^{n+1}`
   is computed from the particles at :math:`{\bf x}_i^{n+1}`.

   Here we can use :math:`{\phi}^n` as an initial guess for :math:`{\phi}^{n+1}` in order to reduce the time
   spent in multigrid to reach the specified tolerance.

Step 7:
   Correct :math:`{\bf U}` with time-centered source terms, and replace :math:`e` by
   :math:`E - \frac{1}{2}U^2` as appropriate.

   We time-center the
   gravitational source terms only,

   .. math:: {\bf U}^{n+1} = {\bf U}^{n+1} + \frac{\Delta t}{2} (S_g^{n+1} - S_g^n)

Step 8:
   Interpolate :math:`{\bf g}^{n+1}` from the grid to the particle locations, then
   update the particle velocities, :math:`{\bf u}_i^{n+1}`

   .. math::

      \begin{aligned}
          {\bf u}_i^{n+1} &=& \frac{1}{a^{n+1}}
                          \left( \left( a^{{n+\frac{1}{2}}} {\bf u}^{{n+\frac{1}{2}}}_i \right)
                               + \frac{\Delta t}{2} \; {\bf g}^{n+1}_i \right)  \end{aligned}

Step \**:
   in post\_ timestep, do a reset and compute\_ new\_ temp after syncing the gravity sources

Deferred-Correction Splitting Algorithm
---------------------------------------

This algorithm is based on the version of SDC introduced by Nonaka et. al. ``\cite{Nonaka2012}``

.. math::

   \frac{\partial \rho}{\partial t} = - \frac{1}{a} \nabla \cdot (\rho {\bf U}) , \\

.. math::

   \begin{aligned}
   \frac{\partial (a \rho {\bf U})}{\partial t} &=& 
   -             \nabla \cdot (\rho {\bf U}{\bf U}) 
   -             \nabla p 
   +             \rho {\bf g} , \end{aligned}

.. math::

   \begin{aligned}
   \frac{\partial (a^2 \rho E)}{\partial t} &=&  a \left[
    -\nabla \cdot (\rho {\bf U}E + p {\bf U})
   +  \rho {\bf U}\cdot {\bf g}
   +  \rho \Lambda_{HC}  \right] . \end{aligned}

.. math::

   \begin{aligned}
   \frac{\partial (a^2 \rho e)}{\partial t} &=& a \left[ 
   - \nabla \cdot (\rho {\bf U}e)
   -  p \nabla \cdot {\bf U}
   +  \rho \Lambda_{HC}  \right] . \end{aligned}

The algorithm at a single level of refinement begins by computing the time step, :math:`\Delta t,`
and advancing :math:`a` from :math:`t^n` to :math:`t^{n+1} = t^n + \Delta t`. The rest of the time step is
then composed of the following steps:

Step 1:
   Compute :math:`{\phi}^n` and :math:`{\bf g}^n` using :math:`\rho^n` and :math:`\rho_{dm}^n`,
   where :math:`\rho_{dm}^{n}` is computed from the particles at :math:`{\bf x}_i^{n}`.

   We note that in the single-level algorithm we can instead use :math:`{\bf g}` as computed at the
   end of the previous step because there have been no changes to :math:`{\bf x}_i`
   or :math:`\rho` since then.

Step 2:
   Interpolate :math:`{\bf g}^n` from the grid to the particle locations, then
   advance the particle velocities by :math:`\Delta t/ 2` and particle positions by :math:`\Delta t`.

   .. math::

      \begin{aligned}
           {\bf u}_i^{{n+\frac{1}{2}}} &=& \frac{1}{a^{{n+\frac{1}{2}}}} ((a^n {\bf u}^n_i) + \frac{\Delta t}{2} \; {\bf g}^n_i) \\
           {\bf x}_i^{n+1}  &=& {\bf x}^n_i + \frac{\Delta t}{a^{{n+\frac{1}{2}}}} {\bf u}_i^{{n+\frac{1}{2}}}\end{aligned}

Step 3:
   Construct advective update terms using :math:`I_R` from last timestep as source

   .. math::

      \begin{aligned}
      A_{\rho} & = & -\frac{1}{a}\nabla\cdot(\rho{\bf U})\\
      A_{\rho u} & = & -\nabla\cdot\left(\rho uu\right)-\nabla p\\
      A_{\rho E} & = & a\left[-\nabla\cdot(\rho{\bf U}E+p{\bf U})\right]\\
      A_{\rho e} & = & \frac{1}{a} \left[
      - \nabla \cdot (\rho_b {\bf U}e)
      - p \nabla \cdot {\bf U}) \right]\end{aligned}

Step 4:
   Update momentum and :math:`\rho E`

   .. math::

      \begin{aligned}
      S_{g} & = & \rho g\\
            &  & \rho{\bf U}\cdot{\bf g}\end{aligned}

   .. math:: u^{n+1,\ast} = u^{n} + \Delta tA^{n+\frac{1}{2}}+ \Delta tS_g^n

   .. math:: \left(\rho E\right)^{n+1,\ast }=\left(\rho E\right)^{n}+ \Delta tA_{\rho E}^{n+1/2} + \Delta tS_g

   After adding the advective update terms :math:`A_\ast` to the appropriate components of :math:`{\bf U}^{n+1,\ast}` and before calculating :math:`S_{g}`, apply a correction to :math:`{\bf U}^{n+1,\ast}` and `A_{\rho}` to enforce :math:`\rho > 1.1 \times small_dens`

Step 5:
   Simultaneously solve heating-cooling:

   .. math::

      \begin{aligned}
      \rho^{n+1,\ast} & = & \rho^{n}+\int_{t^{n}}^{t^{n+1}}A_{\rho}dt^{\prime}\\
      e^{n+1,\ast} & = & e^{n}+\int_{t^{n}}^{t^{n+1}} \left(A_{e}+\Lambda_{HC}\right) dt^{\prime}\end{aligned}

   where :math:`A_{e}=\frac{1}{\Delta t}\left(\left(\left[\frac{1}{a^{n+1}}\right]^{2}\left(\left[a^{n}\right]^{2}\left(\rho e\right)^{n}+\Delta t*A_{\rho e}\right)+A_{reset}\right)/\left(\rho^{n}+\Delta tA_{\rho}\right)-e^{n}\right)`

Step 6:
   We define

   .. math::

      \begin{aligned}
      I_{R_{\left(\rho e\right)}} & = & \left( \left[a^{n+1}\right]^{2}\rho^{n+1,\ast}e^{n+1,\ast}-\left(\left[a^{n}\right]^{2}\rho^{n}e^{n}+\Delta tA_{\rho e}\right)\right)/\left[\Delta t\left(\frac{a^{n}+a^{n+1}}{2}\right)\right]\\
      & & -\left[a^{n+1}\right]^{2}A_{reset}/\left[\Delta t\left(\frac{a^{n}+a^{n+1}}{2}\right)\right]\end{aligned}

Step 7:
   We update internal and total energy using this forcing:
   :math:`\left(\rho e\right)^{n+1,\ast}=\left(\rho e\right)^{n+1,\ast} + \left(\frac{a^{n}+a^{n+1}}{2}\right) \left(\frac{1}{a^{n+1}}\right)^2 \Delta tI_{R_{\rho e}}`
	  
   :math:`\left(\rho E\right)^{n+1,\ast}=\left(\rho E\right)^{n+1,\ast} + \left(\frac{a^{n}+a^{n+1}}{2}\right) \left(\frac{1}{a^{n+1}}\right)^2 \Delta tI_{R_{\rho e}}`

   We store Ne and Temp based on eos\_ hc updates from :math:`(e)^{n+1}`

Step 8:
   Repeat step 3-7

Step 9:
   Compute :math:`{\phi}^{n+1}` and :math:`{\bf g}^{n+1}` using
   :math:`\rho^{n+1,*}` and :math:`\rho_{dm}^{n+1}`, where :math:`\rho_{dm}^{n+1}`
   is computed from the particles at :math:`{\bf x}_i^{n+1}`.

   Here we can use :math:`{\phi}^n` as an initial guess for :math:`{\phi}^{n+1}` in order to reduce the time
   spent in multigrid to reach the specified tolerance.

Step 10:
   Correct :math:`{\bf U}` with time-centered source terms, and replace :math:`e` by
   :math:`E - \frac{1}{2}U^2` as appropriate.

   We time-center the
   gravitational source terms only,

   .. math:: {\bf U}^{n+1} = {\bf U}^{n+1} + \frac{\Delta t}{2} (S_g^{n+1} - S_g^n)

Step 11:
   Interpolate :math:`{\bf g}^{n+1}` from the grid to the particle locations, then
   update the particle velocities, :math:`{\bf u}_i^{n+1}`

   .. math::

      \begin{aligned}
          {\bf u}_i^{n+1} &=& \frac{1}{a^{n+1}}
                          \left( \left( a^{{n+\frac{1}{2}}} {\bf u}^{{n+\frac{1}{2}}}_i \right)
                               + \frac{\Delta t}{2} \; {\bf g}^{n+1}_i \right)  \end{aligned}

Step \**:
   in post\_ timestep, do a reset and compute\_ new\_ temp after syncing the gravity sources
