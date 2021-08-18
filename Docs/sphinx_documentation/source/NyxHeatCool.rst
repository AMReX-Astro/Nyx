*****************************
Radiative Heating and Cooling
*****************************

Nyx provides the capability to compute local heating and cooling effects due to radiation.
The motivation and algorithm for the heating and cooling components is documented in :raw-latex:`\cite{lukic15}`, and the relevant code is located in the ``Source/HeatCool`` subdirectory.
The code is activated through the ``USE_HEATCOOL=TRUE`` option in the ``GNUmakefile``.
Mathematically, the heating and cooling can be described by a single ODE in each cell, to be integrated per time step :math:`\Delta t`.
This ODE exhibits a sensitive relationship to quantities such as temperature and free electron density, and consequently it often requires sophisticated integration techniques to compute correctly.

Nyx provides a few different parallelization strategies techniques for solving this ODE.

For legacy reasons, the integration scheme is selected via the ``nyx.heat_cool_type`` input parameter.
One method is to use ``USE_HEATCOOL=FALSE`` with no ODE solve (selected with ``nyx.heat_cool_type=0``).
The other method is to use ``USE_HEATCOOL=TRUE`` with a vectorized CVODE solve (selected with ``nyx.heat_cool_type=11``).

Users should note that, when compiling with GNUmake , CVODE must be compiled as a separate library; instructions for compiling CVODE are provided in Getting Started.
To link the external CVODE solver into Nyx, one must set ``USE_HEATCOOL=TRUE`` as well as ``USE_SUNDIALS=TRUE`` in the ``GNUmakefile``.

This vectorized option uses CVODE while treating groups of ODEs in different cells as a single system of coupled ODEs.
The purpose of this approach is to enable the evaluation of multiple RHSs simultaneously, for improved parallel performance.
This approach can lead to a significant performance gain in the ODE integration (which is the among the most expensive computational kernels in Nyx).
The group of cells integrated together is controlled by the tilesize used in the CVODE integration loop. If tiling is not used, the entire box of cells is integrated together.
CVODE integrates this resulting system of ODEs with using adaptive time-step control which selects
timesteps for the whole system. We set the absolute tolerance required for the ODE integration in CVODE to ``1e-4``, scaled by the initial value of each cells independent variable in the ODE and
the relative tolerance required for the ODE integration in CVODE to ``1e-4``.
These tolerances, in particular the relative tolerance, have different effects depending on whether one is integrating a single ODE at a time, or a system of ODEs simultaneously.
One should be mindful of the numerical differences which arise from these, which can be observed with the ``fcompare`` tool in AMReX.

Input flags with affect the CVODE integration include:

- ``nyx.use_sundials_constraint`` which when non-zero requires that the internal energy (while the problem is evolved) is positive
- ``nyx.use_sundials_fused`` which when non-zero uses Sundials's GPU fused operations (which are mathematically equivalent, but reduces GPU kernel launch time overhead)
- ``nyx.sundials_alloc_type`` which has up to 5 different vector memory allocation strategies and only affects executables built for GPUs
- ``nyx.use_typical_steps`` which when non-zero sets CVODE's adaptive step-size selection (which substeps the total CVODE integration time) to be ``dt / old_max_steps``. This maximum was over the entire problem domain. (In the strang case, the old maximum is taken from the same phase in the previous time-step)
- ``nyx.sundials_use_tiling`` which controls whether the MFIter loop that iterates over the CVODE integration uses tiling
