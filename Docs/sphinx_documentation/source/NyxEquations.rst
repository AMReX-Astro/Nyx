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
Additionally, :math:`{\bf S}_{\rho {E}}` may include any external forcing terms on the total energy, for example as in the stochastic forcing application.
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

