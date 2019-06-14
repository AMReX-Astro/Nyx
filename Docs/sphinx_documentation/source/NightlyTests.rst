.. _Chap:NightlyTesting :

Nightly Regression Tests
========================

The following regression tests are run nightly with Nyx.   The plotfiles generated in each night's test 
are compared with the benchmark plotfiles using the AMReX :cpp:`fcompare` utility to compare the mesh data
and :cpp:`particle_compare` to compare the particle data.

The results of these tests can be found at https://ccse.lbl.gov/pub/RegressionTesting/Nyx

Below Ng = number of grids, Npa = number of particles, Np = number of MPI ranks / number of OMP threads per rank.
(If just a single number then pure MPI.)

+----------------------+----+----+-------+----+-----+----------------------+
| Test                 | nx | Nl | Npa   | Ng | Np  | What does this test? |
|                      | ny |    |       |    |     |                      |
|                      | nz |    |       |    |     |                      |
+======================+====++====+=======+====+====+======================+
| AMR-density          | 32 | 3  | 5005  | 1  | 4   | Multilevel AMR       |
+----------------------+----+----+-------+----+-----+----------------------+
| DR_restart           | 64 | 3  | 40040 | 8  | 2/2 | Hydro only - restart |
+----------------------+----+----+-------+----+-----+----------------------+
| DoubleRarefaction    | 32 | 3  | 5005  | 8  | 2/2 | Hydro only           |
+----------------------+----+----+-------+----+-----+----------------------+
| DrivenTurbulence     | 10 | 1  | 1611  | 1  | 1/4 | Turbulent forcing    |
+----------------------+----+----+-------+----+-----+----------------------+
| DrivenTurbulence-OMP | 10 | 1  | 1611  | 1  | 1/4 | Forcing with OMP     |
+----------------------+----+----+-------+----+-----+----------------------+
| LyA                  | 10 | 1  | 1611  | 1  | 4   | Includes h/c         |
+----------------------+----+----+-------+----+-----+----------------------+
| Lya-OMP              | 4  | 1  | 2500  | 1  | 1   | Mixed MI/PO + Per    |
+----------------------+----+----+-------+----+-----+----------------------+
| MiniSB               | 4  | 1  | 2500  | 1  | 2   | Mixed MI/PO + Per    |
+----------------------+----+----+-------+----+-----+----------------------+
| MiniSB-ref           | 4  | 2  | 2500  | 1  | 2   | Mixed MI/PO + Per    |
+----------------------+----+----+-------+----+-----+----------------------+
| Part-2line           | 5  | 3  | 1     | 10 | 2   | Single particle      |
+----------------------+----+----+-------+----+-----+----------------------+
| Part-2line_restart   | 5  | 3  | 1     | 10 | 2   | Single particle      |
+----------------------+----+----+-------+----+-----+----------------------+
| Part-mass.nosub      | 5  | 3  | 1     | 10 | 2   | Single particle      |
+----------------------+----+----+-------+----+-----+----------------------+
| Part-mass.sub        | 5  | 3  | 1     | 10 | 2   | Single particle      |
+----------------------+----+----+-------+----+-----+----------------------+
| SantaBarbara         | 5  | 3  | 1     | 10 | 2   | Single particle      |
+----------------------+----+----+-------+----+-----+----------------------+
| Sedov                | 5  | 1  | 1     | 10 | 1   | Hydro only           |
+----------------------+----+----+-------+----+-----+----------------------+
| Sod                  | 5  | 3  | 1     | 10 | 1   | Hydro only           |
+----------------------+----+----+-------+----+-----+----------------------+
| StrongShockTurb      | 5  | 3  | 1     | 10 | 1   | Hydro only           |
+----------------------+----+----+-------+----+-----+----------------------+
| TurbForce            | 5  | 1  | 1     | 10 | 4   | Turbulent forcing    |
+----------------------+----+----+-------+----+-----+----------------------+
| TurbForce-OMP        | 5  | 1  | 1     | 10 | 1   | Forcing with OMP     |
+----------------------+----+----+-------+----+-----+----------------------+

