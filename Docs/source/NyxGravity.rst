*******
Gravity
*******

In Nyx we always compute gravity by solving a Poisson equation on the mesh hierarchy.
 To make sure this option is chosen correctly, we must always set::

   USE_GRAV=TRUE
   
in the GNUmakefile and ::
  
  nyx.do_grav= 1
  
in the inputs file.
To define the gravitational vector we set

  .. math:: \mathbf{g}(\mathbf{x},t) = -\nabla \phi

  where

  .. math:: \mathbf{\Delta} \phi = \frac{4 \pi G}{a} (\rho - \overline{\rho}) \label{eq:Self Gravity}

where :math:`\overline{\rho}` is the average of :math:`\rho` over the entire domain if we assume triply periodic boundary conditions,
and :math:`a(t)` is the scale of the universe as a function of time.
