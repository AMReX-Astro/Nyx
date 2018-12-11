*******
Preface
*******

Welcome to Nyx's User’s Guide!

In this User’s Guide we describe how to download and run Nyx, a massively parallel code
that couples the compressible hydrodynamic equations on a grid with a particle represenation
of dark matter.
Nyx uses an Eulerian grid for the hydrodynamics solver and incorporates adaptive mesh refinement (AMR).Our approach to AMR uses a nested hierarchy of logically-rectangular grids with simultaneous
refinement in both space and time, utilizing the
AMReX library.

For the first paper describing the Nbody, adiabatic hydro and gravity algorithms:
  * A. S. Almgren, J. B. Bell, M.J. Lijewski, Z. Lukic, E. Van Andel, **"Nyx: A Massively Parallel AMR Code for Computational Cosmology"** Astrophysical Journal, **765**, 39, 2013. [[pdf]](../../Publications/almgren/nyx.pdf)
For the reaction and thermal rates of the primordial chemical composition gas (and convergence tests in the context of the Lymanalpha forest), see:
  * Zarija Lukic, Casey Stark, Peter Nugent, Martin White, Avery Meiksin, Ann Almgren, **"The Lymanalpha forest in opticallythin hydrodynamical simulations,"** Monthly Notices of the Royal Astronomical Society, **446**, 36973724, 2015 [[arxiv]](http://arxiv.org/abs/1406.6361).
For the synthesis model of the UV background, which provides the photoionization and photoheating rates, see:
  * Jose Onorbe, Joseph F Hennawi, Zarija Lukic, **"SelfConsistent Modeling of Reionization in Cosmological Hydrodynamical Simulations,"** Astrophysical Journal, **837**, 106, 2017. [[arxiv]](http://arxiv.org/abs/1607.04218)

For Nyx-related Papers:
  * W. Schmidt, C. Byrohl, J. F. Engels, C. Behrens, J. C. Niemeyer, **Viscosity, pressure and support of the gas in simulations of merging coolcore clusters** Monthly Notices of the Royal Astronomical Society, **470(1)**, pp. 142156, 2017 [[arxiv]](http://arxiv.org/abs/1705.06933).  
  * Daniele Sorini, Jose Onorbe, Zarija Lukic, Joseph F Hennawi, **"Modeling the Lymanalpha Forest in Collisionless Simulations,"** Astrophysical Journal, **827**, 97, 2016. [[arxiv]](http://arxiv.org/abs/1602.08099).
  * Brian Friesen, Ann Almgren, Zarija Lukic, Gunther Weber, Dmitriy Morozov, Vincent Beckner, Marcus Day, **"In situ and intransit analysis of cosmological simulations,"**, Computational Astrophysics and Cosmology, **3**, 4, 2016 [[arxiv]](http://arxiv.org/abs/1705.06933).
  * W. Schmidt, J.F. Engels, J.C. Niemeyer, A.S. Almgren, **"Hot and Turbulent Gas in Clusters",** Monthly Notices of the Royal Astronomical Society, **459(1)**, pp. 701719, 2016 [[arxiv]](http://arxiv.org/abs/1603.04711).
  * H. Braun, W. Schmidt, J.C. Niemeyer, A.S. Almgren, **"Largeeddy simulations of isolated disk galaxies with thermal and turbulent feedback,"** Monthly Notices of the Royal Astronomical Society, **442(4)**, pp. 34073426, 2014. [[arxiv]](http://arxiv.org/pdf/1405.6245)
  * W. Schmidt, A.S. Almgren, H. Braun, J.F. Engels, J.C. Niemeyer, R.R. Mekuria, A.J. Aspden, J.B. Bell, **"Cosmological Fluid Mechanics with Adaptively Refined Large Eddy Simulations,"** Monthly Notices of the Royal Astronomical Society, **440**, pp. 30513077, 2014. [[arxiv]](http://arxiv.org/pdf/1309.3996v1.pdf)
  * W. Schmidt, J. Schulz, L. Iapichino, A.S. Almgren, **" Influence of adaptive mesh refinement and the hydro solver on shearinduced mass stripping in a minor merger scenario,"** Astronomy and Computing, **9**, pp. 4964, March 2015 [[arxiv]](http://arxiv.org/abs/1411.7275).
		  
The development of AMReX library is led by the
Center for Computational Sciences and Engineering / Lawrence Berkeley
National Laboratory. Nyx development is done collaboratively,
including the Computational Cosmology Center and CCSE. All of Nyx's development is done in the public github repository—anyone can see the updates as they are done.  We are always happy to have new developers as part of the Nyx team. Fork the Nyx git repository on github, make your changes, and issue a pull request against the master github repo. Any level of changes are welcomed: documentation, bug fixes, new test problems, new solvers, ...

 .. todo::
    Describe developer/contributor/author list

