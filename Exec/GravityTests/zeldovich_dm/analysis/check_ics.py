"""
Script for analyzing the Zeldovich initial conditions. Expects files
``zeldovich_particles.ascii`` and ``zeldovich_params.csv``.

Verifies that the normal vector is a unit vector and that the velocities along
the sheets are 0.

Plots the particles in a 3D projection.

Author: Casey W. Stark <caseywstark@gmail.com>
Copyright 2011.
Included September 19, 2011.

"""

import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

from zeldovich import get_scale_factor, perturbed_x, perturbed_v_x

grav_const = 6.673e-8;  # in cm^3 g^-1 s^-2
# units factor for critical density (rho_c), km^2 Mpc^-2 cm^-3 g -> M_sun Mpc^-3
# terms are (km^2 / cm^2) * (Mpc / cm) * (g / M_sun)
rho_c_units_factor = 1.0e10 * 3.08568025e24 * 5.02785431e-34

pass_count = 0
fail_count = 0
num_tests = 4
test_index = 0

print ""
print "Testing Zeldovich initial conditions."

# path to initial conditions files
ic_path = os.path.abspath("../ics/")

### Read in params
reader = csv.DictReader(open(os.path.join(ic_path, "zeldovich_params.csv"), "r"))
params = reader.next()

initial_z = float(params["initial_z"])
caustic_z = float(params["caustic_z"])
current_a = get_scale_factor(initial_z)
box_length = float(params["box_length"])
k = float(params["k"])
offset = float(params["offset"])
H_0 = float(params["H_0"])
sheet_normal = np.array(map(float, (params["normal_x"], params["normal_y"],
                                    params["normal_z"])))

### Read in particle data
sys.stdout.write("Reading in particle data ... ")
sys.stdout.flush()

p_data = np.loadtxt(os.path.join(ic_path, "zeldovich_particles.ascii"), skiprows=1)

print "DONE"

sys.stdout.write("Processing particle data ... ")
sys.stdout.flush()

# reference columns into aliases and rotate.
# note we divide the velocities by the scale factor to make them comoving.
xs = p_data[:, 0]
ys = p_data[:, 1]
zs = p_data[:, 2]
masses = p_data[:, 3]
num_particles = len(xs)
# Distances from the plane defined by the normal at (0, 0, 0)
pos_normals = (xs * sheet_normal[0] + ys * sheet_normal[1]
               + zs * sheet_normal[2])

# velocities, make them comoving
vxs = p_data[:, 4] / current_a
vys = p_data[:, 5] / current_a
vzs = p_data[:, 6] / current_a

# Get the velocities out of the plane and in the plane.
v_normals = (vxs * sheet_normal[0] + vys * sheet_normal[1]
             + vzs * sheet_normal[2])

v_alongs = np.sqrt(np.power(vys * sheet_normal[2] - vzs * sheet_normal[1], 2)
                   + np.power(vzs * sheet_normal[0] - vxs * sheet_normal[2], 2)
                   + np.power(vxs * sheet_normal[1] - vys * sheet_normal[0], 2))

print "DONE"
print ""
print "Running tests..."

test_index += 1

### Check that sheet_normal is a unit vector.
sys.stdout.write("(%i/%i) Is sheet normal vector a unit vector ... "
                 % (test_index, num_tests))
sys.stdout.flush()
norm = np.linalg.norm(sheet_normal)
if norm == 1.0:
    print "PASS"
    pass_count += 1
else:
    print "FAIL"
    print "Normal vector is %s" % sheet_normal
    print "with norm %f" % norm
    fail_count += 1

test_index += 1

### Check v along is 0
sys.stdout.write("(%i/%i) Is velocity along sheets 0 ... "
                 % (test_index, num_tests))
sys.stdout.flush()
good_v_alongs = (v_alongs == 0.0)
good_v_along = good_v_alongs.all()
if good_v_along:
    print "PASS"
    pass_count += 1
else:
    print "FAIL"
    fail_count += 1

test_index += 1

### Check that box mass density is consistent with cosmology
density_tolerance = 1  # M_sun / Mpc^3
sys.stdout.write("(%i/%i) H_0 value gives critical density matching density ...\n"
                 % (test_index, num_tests))
sys.stdout.flush()
critical_density = (3 * H_0**2 / (8 * np.pi * grav_const) * rho_c_units_factor)
box_density = masses.sum() / box_length**3
density_residual = abs(box_density - critical_density)

print ""
print "Critical density:\t%.4e" % critical_density
print "Box density:\t\t%.4e" % box_density
print "Residual:\t\t%.4e" % density_residual

if density_residual <= density_tolerance:
    print "PASS"
    pass_count += 1
else:
    print "FAIL"
    fail_count += 1
print ""

test_index += 1

### Check that none of the particles sit on the boundary
sys.stdout.write("(%i/%i) No particles sit on boundary or outside box ... "
                 % (test_index, num_tests))
sys.stdout.flush()

bad_xs_0 = (xs <= 0.0).any()
bad_xs_l = (xs >= box_length).any()
bad_ys_0 = (ys <= 0.0).any()
bad_ys_l = (ys >= box_length).any()
bad_zs_0 = (zs <= 0.0).any()
bad_zs_l = (zs >= box_length).any()

bad_pos = False
if bad_xs_0:
    print "Some x coordinate at 0 or smaller!"
    bad_pos = True
if bad_xs_l:
    print "Some x coordinate at L or larger!"
    bad_pos = True
if bad_ys_0:
    print "Some y coordinate at 0 or smaller!"
    bad_pos = True
if bad_ys_l:
    print "Some y coordinate at L or larger!"
    bad_pos = True
if bad_zs_0:
    print "Some z coordinate at 0 or smaller!"
    bad_pos = True
if bad_zs_l:
    print "Some z coordinate at L or larger!"
    bad_pos = True

if not bad_pos:
    print "PASS"
    pass_count += 1
else:
    print "FAIL"
    fail_count += 1

test_index += 1

# Report
print ""
print "Passed: %i    Failed: %i" % (pass_count, fail_count)

### Check if we match Zeldovich predictions
unique_ds = np.unique(pos_normals)
normal_length = min(box_length / sheet_normal[0],
                    box_length / sheet_normal[1],
                    box_length / sheet_normal[2])

q_array = np.linspace(0, normal_length, 100 * normal_length * k)

# Zeldovich predictions
zel_xs = perturbed_x(q_array + offset, initial_z, caustic_z, k) - offset
zel_vxs = perturbed_v_x(q_array + offset, initial_z, caustic_z, k, H_0)


### Visual check: 3D matplotlib projection of particle positions
print ""
sys.stdout.write("Making plots ... ")
sys.stdout.flush()

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_title("Zeldovich phase space")
ax.set_xlabel("Comoving normal position (Mpc)")
ax.set_ylabel("Proper normal velocity (km/s)")
ax.plot(pos_normals, v_normals, color="None", linestyle="None", marker=".",
        markerfacecolor=None, markeredgecolor="red", alpha=0.2)
ax.plot(zel_xs, zel_vxs, color="gray", linestyle="-")

plot_name_1 = "ICs_phase_space.png"
fig.savefig(plot_name_1)
plt.close()

# particle positions
arrow_size = normal_length  #box_length / 10.0
arrow_x = np.array((0, arrow_size * sheet_normal[0]))
arrow_y = np.array((0, arrow_size * sheet_normal[1]))
arrow_z = np.array((0, arrow_size * sheet_normal[2]))

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_title("Zeldovich Particle Positions")
ax.set_xlabel("x (Mpc)")
ax.set_ylabel("y (Mpc)")
ax.set_zlabel("z (Mpc)")
ax.view_init(10, -80)
ax.plot(arrow_x, arrow_y, arrow_z, '-')
ax.plot(xs, ys, zs, '.', alpha=0.02)

plot_name_2 = "ICs_3D.png"
fig.savefig(plot_name_2)
plt.close()

print "DONE"
print "Saved plots as `%s` and `%s`" % (plot_name_1, plot_name_2)
print ""
