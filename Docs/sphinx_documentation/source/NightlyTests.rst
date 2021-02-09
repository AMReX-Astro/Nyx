.. _Chap:NightlyTesting :

Nightly Regression Tests
========================

The following regression tests are run nightly with Nyx.   The plotfiles generated in each night's test 
are compared with the benchmark plotfiles using the AMReX :cpp:`fcompare` utility to compare the mesh data
and :cpp:`particle_compare` to compare the particle data.

The results of these tests can be found at https://ccse.lbl.gov/pub/RegressionTesting/Nyx.
These tests are also run on an NVIDIA GPU and those results can be found at https://ccse.lbl.gov/pub/GpuRegressionTesting/Nyx.

We have a number of compile-time options -- this chart summarizes which regression tests
test which compile-time options.  (Note that the USE_GRAV and AMREX_USE_PARTICLES macros
which are used in the code source are both set to TRUE in the build system.)

+---------------------------+----------+--------------+
| Test                      | NO_HYDRO | USE_HEATCOOL |
+===========================+==========+==============+
| AMR-density               | FALSE    |   TRUE       |
+---------------------------+----------+--------------+
| DR_restart                | FALSE    |   FALSE      |
+---------------------------+----------+--------------+
| DoubleRarefaction         | FALSE    |   FALSE      |
+---------------------------+----------+--------------+
| DrivenTurbulence          | FALSE    |   FALSE      |
+---------------------------+----------+--------------+
| DrivenTurbulence-OMP      | FALSE    |   FALSE      |
|  (CPU only)               |          |              |
+---------------------------+----------+--------------+
| LyA                       | FALSE    |   TRUE       |
+---------------------------+----------+--------------+
| LyA-adiabatic             | FALSE    |   FALSE      |
+---------------------------+----------+--------------+
| LyA-OMP                   | FALSE    |   TRUE       |
|  (CPU only)               |          |              |
+---------------------------+----------+--------------+
| LyA_Neutrinos             | FALSE    |   TRUE       |
|  (CPU only)               |          |              |
+---------------------------+----------+--------------+
| MiniSB                    | FALSE    |  FALSE       |
+---------------------------+----------+--------------+
| MiniSB-ppm                | FALSE    |  FALSE       |
+---------------------------+----------+--------------+
| MiniSB-ref                | FALSE    |  FALSE       |
+---------------------------+----------+--------------+
| MiniSB-ref-ppm            | FALSE    |  FALSE       |
+---------------------------+----------+--------------+
| Part-2line                | FALSE    |  FALSE       |
+---------------------------+----------+--------------+
| Part-2line-nbody          | FALSE    |  FALSE       |
+---------------------------+----------+--------------+
| Part-2line-nomesh         | TRUE     |  FALSE       |
+---------------------------+----------+--------------+
| Part-2line_restart        | FALSE    |  FALSE       |
+---------------------------+----------+--------------+
| Part-2line_restart-nbody  | FALSE    |  FALSE       |
+---------------------------+----------+--------------+
| Part-2line_restart-nomesh | TRUE     |  FALSE       |
+---------------------------+----------+--------------+
| Part-mass.nosub           | FALSE    |  FALSE       |
+---------------------------+----------+--------------+
| Part-mass.nosub-nbody     | FALSE    |  FALSE       |
+---------------------------+----------+--------------+
| Part-mass.nosub-nomesh    | TRUE     |  FALSE       |
+---------------------------+----------+--------------+
| Part-mass.sub             | FALSE    |  FALSE       |
+---------------------------+----------+--------------+
| Part-mass.sub-nbody       | FALSE    |  FALSE       |
+---------------------------+----------+--------------+
| Part-mass.sub-nbody       | FALSE    |  FALSE       |
+---------------------------+----------+--------------+
| Part-mass.sub-nomesh      | TRUE     |  FALSE       |
+---------------------------+----------+--------------+
| Sedov                     | FALSE    |  FALSE       |
+---------------------------+----------+--------------+
| Sedov-ppm                 | FALSE    |  FALSE       |
+---------------------------+----------+--------------+
| Sod                       | FALSE    |  FALSE       |
+---------------------------+----------+--------------+
| Sod-ppm                   | FALSE    |  FALSE       |
+---------------------------+----------+--------------+
| StrongShockTube           | FALSE    |  FALSE       |
+---------------------------+----------+--------------+
| StrongShockTube-ppm       | FALSE    |  FALSE       |
+---------------------------+----------+--------------+

Below nx,ny,nz are the number of cells in each coordinate direction at the coarsest level;
Npa = number of particles, Np = number of MPI ranks / number of OMP threads per rank.
(If just a single number then pure MPI.)

+---------------------------+----------+----+--------+-----+----------------------+
| Test                      | nx ny nz | Nl | Npa    | Np  | What does this test? |
+===========================+==========+====+========+=====+======================+
| AMR-density               | 64 64 64 | 3  | 262144 | 4   | Start from chk00300  |
+---------------------------+----------+----+--------+-----+----------------------+
| DR_restart                | 32  4  4 | 3  | 0      | 2/2 | Hydro only - restart |
+---------------------------+----------+----+--------+-----+----------------------+
| DoubleRarefaction         | 32  4  4 | 3  | 0      | 2/2 | Hydro only           |
+---------------------------+----------+----+--------+-----+----------------------+
| DrivenTurbulence          | 32 32 32 | 1  | 0      | 4   | Turbulent forcing    |
+---------------------------+----------+----+--------+-----+----------------------+
| DrivenTurbulence-OMP      | 32 32 32 | 1  | 0      | 1/4 |  Turbulent forcing   |
|  (CPU only)               |          |    |        |     |  with OMP            |
+---------------------------+----------+----+--------+-----+----------------------+
| LyA                       | 32 32 32 | 1  | 32768  | 4   | Includes h/c         |
+---------------------------+----------+----+--------+-----+----------------------+
| Lya-OMP                   | 32 32 32 | 1  | 32768  | 1/4 | LyA with OMP         |
|  (CPU only)               |          |    |        |     |                      |
+---------------------------+----------+----+--------+-----+----------------------+
| LyA-adiabatic             | 32 32 32 | 1  | 32768  | 4   | No h/c               |
+---------------------------+----------+----+--------+-----+----------------------+
| LyA_Neutrinos             | 32 32 32 | 1  | 32768  | 1   | LyA with OMP         |
|  (CPU only)               |          |    |        |     |                      |  
+---------------------------+----------+----+--------+-----+----------------------+
| MiniSB                    | 32 32 32 | 1  | 2500   | 2   | Small version of SB  |
+---------------------------+----------+----+--------+-----+----------------------+
| MiniSB-ppm                | 32 32 32 | 1  | 2500   | 2   | Small version of SB  |
+---------------------------+----------+----+--------+-----+----------------------+
| MiniSB-ref                | 32 32 32 | 2  | 2500   | 2   | Small version of SB  |
|                           |          |    |        |     | with refinement      |
+---------------------------+----------+----+--------+-----+----------------------+
| MiniSB-ref-ppm            | 32 32 32 | 2  | 2500   | 2   | Small version of SB  |
|                           |          |    |        |     | with refinement      |
+---------------------------+----------+----+--------+-----+----------------------+
| Part-2line                | 16 16 16 | 3  | 16     | 2/2 | Particle-only        |
+---------------------------+----------+----+--------+-----+----------------------+
| Part-2line-nbody          | 16 16 16 | 3  | 16     | 2/2 | Particle-only        |
+---------------------------+----------+----+--------+-----+----------------------+
| Part-2line-nomesh         | 16 16 16 | 3  | 16     | 2/2 | Particle-only        |
+---------------------------+----------+----+--------+-----+----------------------+
| Part-2line_restart        | 16 16 16 | 3  | 16     | 2/2 | Particle-only        |
+---------------------------+----------+----+--------+-----+----------------------+
| Part-2line_restart-nbody  | 16 16 16 | 3  | 16     | 2/2 | Particle-only        |
+---------------------------+----------+----+--------+-----+----------------------+
| Part-2line_restart-nomesh | 16 16 16 | 3  | 16     | 2/2 | Particle-only        |
+---------------------------+----------+----+--------+-----+----------------------+
| Part-mass.nosub           | 16 16 16 | 3  | 8      | 2/2 | Particle-only        |
+---------------------------+----------+----+--------+-----+----------------------+
| Part-mass.nosub-nbody     | 16 16 16 | 3  | 8      | 2/2 | Particle-only        |
+---------------------------+----------+----+--------+-----+----------------------+
| Part-mass.nosub-nomesh    | 16 16 16 | 3  | 8      | 2/2 | Particle-only        |
+---------------------------+----------+----+--------+-----+----------------------+
| Part-mass.sub             | 16 16 16 | 3  | 8      | 2/2 | Particle-only        |
+---------------------------+----------+----+--------+-----+----------------------+
| Part-mass.sub-nbody       | 16 16 16 | 3  | 8      | 2/2 | Particle-only        |
+---------------------------+----------+----+--------+-----+----------------------+
| Part-mass.sub-nomesh      | 16 16 16 | 3  | 8      | 2/2 | Particle-only        |
+---------------------------+----------+----+--------+-----+----------------------+
| Sedov                     | 32 32 32 | 1  | 0      | 1   | Hydro only (PLM)     |
+---------------------------+----------+----+--------+-----+----------------------+
| Sedov-ppm                 | 32 32 32 | 1  | 0      | 1   | Hydro only (PPM)     |
+---------------------------+----------+----+--------+-----+----------------------+
| Sod                       | 4  32 4  | 3  | 0      | 1   | Hydro only (PLM)     |
+---------------------------+----------+----+--------+-----+----------------------+
| Sod-ppm                   | 4  32 4  | 3  | 0      | 1/2 | Hydro only (PPM)     |
+---------------------------+----------+----+--------+-----+----------------------+
| StrongShockTube           | 32 4  4  | 3  | 0      | 1/2 | Hydro only (PLM)     |
+---------------------------+----------+----+--------+-----+----------------------+
| StrongShockTube-ppm       | 32 4  4  | 3  | 0      | 1   | Hydro only (PPM)     |
+---------------------------+----------+----+--------+-----+----------------------+

