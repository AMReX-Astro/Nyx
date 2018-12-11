==============================================
Hydrodynamic Equations in Comoving Coordinates
==============================================

Conservative Form
-----------------

We solve the equations of gas dynamics in a coordinate system that is comoving
with the expanding universe, with expansion factor, :math:`a,` related to the redshift, :math:`z`, by :math:`a = 1 / (1 + z).`
The continuity equation is written,

.. math::

   \label{eq:dens}
   \frac{\partial \rho_b}{\partial t} = - \frac{1}{a} \nabla \cdot (\rho_b {\bf U})  , \\

where :math:`\rho_b` is the comoving baryonic density, related to the proper density by :math:`\rho_b = a^3 \rho_{proper},`
and :math:`{\bf U}` is the proper peculiar baryonic velocity.

The momentum evolution equation can be expressed as

.. math::

   \begin{aligned}
   \frac{\partial (\rho_b {\bf U})}{\partial t} &=&  \frac{1}{a} \left(
   - \nabla \cdot (\rho_b {\bf U} {\bf U}) 
   - \nabla p 
   + \rho_b {\bf g} 
   + {\bf S}_{\rho {\bf U}}
   - \dot{a} \rho_b {\bf U} \right)  , \end{aligned}

or equivalently,

.. math::

   \begin{aligned}
   \label{eq:momt}
   \frac{\partial (a \rho_b {\bf U})}{\partial t} &=& 
   -             \nabla \cdot (\rho_b {\bf U} {\bf U}) 
   -             \nabla p 
   +             \rho_b {\bf g} 
   +             {\bf S}_{\rho {\bf U}}  , \end{aligned}

where the pressure, :math:`p`, that appears in the
evolution equations is related to the proper pressure, :math:`p_{proper},` by :math:`p = a^3 p_{proper}.`
Here :math:`{\bf g} = - \nabla \phi` is the gravitational acceleration vector, and
:math:`{\bf S}_{\rho {\bf U}}` represents any external forcing terms.

The energy equation can be written,

.. math::

   \begin{aligned}
   \frac{\partial (\rho_b E)}{\partial t} &=& \frac{1}{a} \left[
   - \nabla \cdot (\rho_b {\bf U} E + p {\bf U})
   + ( \rho_b {\bf U} \cdot {\bf g} +  S_{\rho E} ) 
   - \dot{a} ( 3 (\gamma - 1) \rho_b e + \rho_b ( {\bf U} \cdot {\bf U}) ) \right]  . \end{aligned}

or equivalently,

.. math::

   \begin{aligned}
   \label{eq:energy}
   \frac{\partial (a^2 \rho_b E)}{\partial t} &=& a \left[
   - \nabla \cdot (\rho_b {\bf U} E + p {\bf U})
   +  \rho_b {\bf U} \cdot {\bf g} 
   +  S_{\rho E}  
   +  \dot{a} ( \; ( 2 - 3 (\gamma - 1) ) \; \rho_b e ) \right]  . \end{aligned}

Here :math:`E = e + {\bf U} \cdot {\bf U} / 2` is the total energy per unit mass,
where :math:`e` is the specific internal energy.
:math:`S_{\rho E} = S_{\rho e} + {\bf U} \cdot {\bf S}_{\rho {\bf U}}`
where :math:`S_{\rho e} = \Lambda^H - \Lambda^C` represents the heating and cooling terms, respectively.
We can write the evolution equation for internal energy as

.. math::

   \begin{aligned}
   \frac{\partial (\rho_b e)}{\partial t} &=& \frac{1}{a} \left[
   - \nabla \cdot (\rho_b {\bf U} e)
   - p \nabla \cdot {\bf U}
   - \dot{a} ( 3 (\gamma - 1) \rho_b e )
   + S_{\rho e}  \right]  . \end{aligned}

or equivalently,

.. math::

   \begin{aligned}
   \frac{\partial (a^2 \rho_b e)}{\partial t} &=&  a \left[
   - \nabla \cdot (\rho_b {\bf U} e)
   - p \nabla \cdot {\bf U}
   + S_{\rho e} 
   + \dot{a} ( \; ( 2 - 3 (\gamma - 1) ) \; \rho_b e ) \right]  . \end{aligned}

Note that for a gamma-law gas with :math:`\gamma = 5/3,` we can write

.. math::

   \begin{aligned}
   \frac{\partial (a^2 \rho_b E)}{\partial t} &=&  a \left[
    -\nabla \cdot (\rho_b {\bf U} E + p {\bf U})
   +  \rho_b {\bf U} \cdot {\bf g} 
   +  S_{\rho e}  \right]   . \end{aligned}

and

.. math::

   \begin{aligned}
   \frac{\partial (a^2 \rho_b e)}{\partial t} &=& a \left[ 
   - \nabla \cdot (\rho_b {\bf U} e)
   -  p \nabla \cdot {\bf U}
   +  S_{\rho e}  \right]   . \end{aligned}

Tracing
-------

In order to compute the fluxes on faces, we trace :math:`\rho, {\bf U}, \rho e` and :math:`p` to the faces.

Thus we must convert the momentum evolution equation into a velocity evolution equation:

.. math::

   \begin{aligned}
   \frac{\partial{\bf U}}{\partial t} &=&  \frac{1}{\rho_b} \left(
   \frac{\partial (\rho_b {\bf U})}{\partial t}  - {\bf U} \frac{\partial \rho_b}{\partial t}  \right) \\
   &=&  \frac{1}{a \rho_b} \left(
   - \nabla \cdot (\rho_b {\bf U} {\bf U})
   - \nabla p
   + \rho_b {\bf g}
   + S_{\rho {\bf U}}
   - \dot{a} \rho_b {\bf U} \right) 
   + \frac{1}{a}  {\bf U} \; \nabla \cdot (\rho_b {\bf U}) \\
   &=&  \frac{1}{a} \left(
   - {\bf U} \cdot \nabla {\bf U}
   - \frac{1}{\rho_b} \nabla p
   + {\bf g}
   + \frac{1}{\rho_b} {\bf S}_{\rho {\bf U}}
   - \dot{a} {\bf U} \right)  .\end{aligned}

Subgrid Scale Model in Comoving Coordinates
===========================================

The fundamental modification to the standard compressible equations is the addition
of a SGS turbulence energy variable, :math:`K` and associated source terms in the equations
for the evolution of velocity, total energy, and :math:`K` :raw-latex:`\cite{schumann:1975,sagaut,schmidt:2006}`.
The set of conservation equations in comoving coordinates (\ `[eq:dens] <#eq:dens>`__)–(\ `[eq:energy] <#eq:energy>`__) then becomes
:raw-latex:`\cite{maier:2009}`:

.. math::

   \begin{aligned}
   \frac{\partial \rho_b}{\partial t} =& - \frac{1}{a} \nabla \cdot (\rho_b {\bf U})  , \\
   \label{eq:momt_les}
   \frac{\partial (a \rho_b {\bf U})}{\partial t} =& 
   -             \nabla \cdot (\rho_b {\bf U} {\bf U}) 
   -             \nabla p
   +             \nabla \cdot \boldsymbol{\tau}
   +             \rho_b {\bf g}   , \\
   \label{eq:energy_les}
   \frac{\partial (a^2 \rho_b E)}{\partial t} =& - a \nabla \cdot (\rho_b {\bf U} E + p {\bf U})
   +  a \rho_b {\bf U} \cdot {\bf g} 
   + a \nabla \cdot ({\bf U}\cdot\boldsymbol{\tau}) - a^2(\Sigma - \rho_b \varepsilon)  \\
   \nonumber
   &+ a \dot{a} \left( ( 2 - 3 (\gamma - 1) ) \rho_b e \right) 
   + a^2 ( \Lambda^H  - \Lambda^C )   , \\
   \label{eq:k_les}
   \frac{\partial (a^2\rho_b K)}{\partial t} =&
   - a\nabla \cdot \left(\rho_b {\bf U} K\right) 
   + a\nabla \cdot \left(\rho_b \kappa_{\rm sgs}\nabla K\right)  
   + a^2(\Sigma - \rho_b \varepsilon) .\end{aligned}

The interaction between resolved and unresolved turbulent eddies is described by
the SGS turbulence stress tensor :math:`\boldsymbol{\tau}`. Since inertial-range dynamics of
turbulence is scale-invariant, we conjecture that :math:`\boldsymbol{\tau}` in comoving coordinates
has the same form as for non-expanding fluids. For compressible turbulence,
the following closure is proposed in :raw-latex:`\cite{schmidt:2011}`:

.. math::

   \label{eq:tau_nonlin}
     \tau_{ij}= 2C_{1}\Delta\rho_b(2 K_{\mathrm{sgs}})^{1/2}S_{\! ij}^{\ast}
     -4C_{2}\rho_bK\frac{U_{i,k}U_{j,k}}{|\nabla{\bf U}|^{2}}
     -\frac{2}{3}(1-C_{2})\rho_bK\delta_{ij}.

where :math:`|\nabla{\bf U}|:=(2U_{i,k}U_{i,k})^{1/2}` is the norm of the resolved velocity derivative,

.. math::

   S_{ij}^{\ast} = S_{ij} - \frac{1}{3}\delta_{ij}d =
   \frac{1}{2} (U_{i,j} + U_{j,i}) - \frac{1}{3}\delta_{ij}U_{k,k}

is the trace-free rate-of strain, and :math:`\Delta=(\vartriangle\!x\,\vartriangle\!y\,\vartriangle\!z)^{1/3}`
is the grid scale in comoving coordinates. The production and dissipation terms in equation (\ `[eq:k_les] <#eq:k_les>`__)
are defined as follows:

.. math::

   \begin{aligned}
   \Sigma             &=& \frac{1}{a}\tau_{ij} S_{ij}, \\
   \varepsilon &=& \frac{C_\varepsilon K^{3/2}}{a\Delta},\end{aligned}

and :math:`\kappa_{\rm sgs} = C_{\kappa}\Delta K^{1/2}` is the SGS diffusivity.
Here we assume that the Reynolds number of turbulence is high such that the damping of turbulent eddies by the microscopic
viscosity of the fluid occurs entirely on the subgrid scales. Because of the numerical viscosity of PPM,
however, part of the numerically resolved kinetic energy will be dissipated directly into internal energy.
