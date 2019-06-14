.. _Chap:NightlyTesting :

Nightly Regression Tests
========================

The following regression tests are run nightly with Nyx.   The plotfiles generated in each night's test 
are compared with the benchmark plotfiles using the AMReX :cpp:`fcompare` utility to compare the mesh data
and :cpp:`particle_compare` to compare the particle data.

The results of these tests can be found at https://ccse.lbl.gov/pub/RegressionTesting/Nyx

Below nx,ny,nz are the number of cells in each coordinate direction at the coarsest level;
Npa = number of particles, Np = number of MPI ranks / number of OMP threads per rank.
(If just a single number then pure MPI.)

+----------------------+----------+----+--------+-----+----------------------+
| Test                 | nx ny nz | Nl | Npa    | Np  | What does this test? |
+======================+==========+====+========+=====+======================+
| AMR-density          | 64 64 64 | 3  | 262144 | 4   | Start from chk00300  |
+----------------------+----------+----+--------+-----+----------------------+
| DR_restart           | 32  4  4 | 3  | 0      | 2/2 | Hydro only - restart |
+----------------------+----------+----+--------+-----+----------------------+
| DoubleRarefaction    | 32  4  4 | 3  | 0      | 2/2 | Hydro only           |
+----------------------+----------+----+--------+-----+----------------------+
| DrivenTurbulence     | 32 32 32 | 1  | 0      | 1/4 | Turbulent forcing    |
+----------------------+----------+----+--------+-----+----------------------+
| DrivenTurbulence-OMP | 32 32 32 | 1  | 0      | 1/4 | Forcing with OMP     |
+----------------------+----------+----+--------+-----+----------------------+
| LyA                  | 32 32 32 | 1  | 32768  | 4   | Includes h/c         |
+----------------------+----------+----+--------+-----+----------------------+
| Lya-OMP              | 32 32 32 | 1  | 32768  | 1   | Mixed MI/PO + Per    |
+----------------------+----------+----+--------+-----+----------------------+
| MiniSB               | 32 32 32 | 1  | 2500   | 2   | Mixed MI/PO + Per    |
+----------------------+----------+----+--------+-----+----------------------+
| MiniSB-ref           | 32 32 32 | 2  | 2500   | 2   | Mixed MI/PO + Per    |
+----------------------+----------+----+--------+-----+----------------------+
| Part-2line           | 16 16 16 | 3  | 16     | 2   | Particle-only        |
+----------------------+----------+----+--------+-----+----------------------+
| Part-2line_restart   | 16 16 16 | 3  | 16     | 2   | Particle-only        |
+----------------------+----------+----+--------+-----+----------------------+
| Part-mass.nosub      | 16 16 16 | 3  | 8      | 2   | Particle-only        |
+----------------------+----------+----+--------+-----+----------------------+
| Part-mass.sub        | 16 16 16 | 3  | 8      | 2   | Particle-only        |
+----------------------+----------+----+--------+-----+----------------------+
| SantaBarbara         |          |    |        |     | Compile-only         |
+----------------------+----------+----+--------+-----+----------------------+
| Sedov                | 32 32 32 | 1  | 0      | 1   | Hydro only           |
+----------------------+----------+----+--------+-----+----------------------+
| Sod                  | 4  32 4  | 3  | 0      | 1   | Hydro only           |
+----------------------+----------+----+--------+-----+----------------------+
| StrongShockTurb      | 32 4  4  | 3  | 0      | 1   | Hydro only           |
+----------------------+----------+----+--------+-----+----------------------+
| TurbForce            | 32 32 32 | 1  | 0      | 4   | Turbulent forcing    |
+----------------------+----------+----+--------+-----+----------------------+
| TurbForce-OMP        | 32 32 32 | 1  | 0      | 1/4 | Forcing with OMP     |
+----------------------+----------+----+---------+-----+----------------------+

