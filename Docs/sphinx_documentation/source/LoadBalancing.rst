.. role:: cpp(code)
   :language: c++

.. _sec:load_balancing:

Load Balancing
--------------

The process of load balancing is typically independent of the process of grid creation; 
the inputs to load balancing are a given set of grids with a set of weights 
assigned to each grid.

Single-level load balancing algorithms are sequentially applied to each AMR level independently, 
and the resulting distributions are mapped onto the ranks taking into account the weights 
already assigned to them (assign heaviest set of grids to the least loaded rank)

Options supported by AMReX include the following; the default is SFC:

- Knapsack: the default weight of a grid in the knapsack algorithm is the number of grid cells, 
  but AMReX supports the option to pass an array of weights – one per grid – or alternatively 
  to pass in a MultiFab of weights per cell which is used to compute the weight per grid

- SFC: enumerate grids with a space-filling Z-morton curve, then partition the 
  resulting ordering across ranks in a way that balances the load

- Round-robin: sort grids and assign them to ranks in round-robin fashion -- specifically
  FAB *i* is owned by CPU *i* %N where N is the total number of MPI ranks.

Load Balancing the Hydrodynamical Mesh
--------------------------------------

For Nyx, the DistributionMapping defaults to using a SFC strategy.
Setting the DistributionMapping to use a different strategy is controlled
by flags beginning with DistributionMapping::

  DistributionMapping.strategy = {KNAPSACK, SFC, ROUNDROBIN}
  DistributionMapping.verbose = {0, 1}

These flags take effect whenever the Regrid operation is called on the mesh,
typically either when ``amr.regrid_int`` is reached in a multilevel simulation or if
``amr.regrid_on_restart=1``.

Load Balancing the Dark Matter Particles
----------------------------------------

For Nyx, the flags which affect how the active particles get their DistributionMapping change are::

  nyx.load_balance_int = -1
  nyx.load_balance_start_z = 15
  nyx.load_balance_wgt_strategy = {0, 1, 2}
  nyx.load_balance_wgt_nmax = -1
  nyx.load_balance_strategy = {KNAPSACK, SFC, ROUNDROBIN}

The ``wgt_strategy`` uses either the 0=number of cells, 1=number of dark matter particles, or 2=total count of dark matter particles to determine the cost of each box. 1 and 2 should be equivalent, and 1 should be cheaper.

This load balancing is not on by default, and has a flag to determine based on redshift z when to begin affecting the calculation because this introduces additional communication penalties when the particles interact with mesh data. This additional communication is a trade-off to balance the amount of particle memory and particle work between ranks.

This strategy is not currently implemented to work for multilevel simulations.
