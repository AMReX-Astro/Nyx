"""
Script for analyzing the Zeldovich test case output. Expects a file
``zeldovich_params.csv`` in the current dir and  ``particles.ascii`` in the
plotfile directories.

Checks:
 - Velocities along the sheet are small.
 - Phase space.

Also does density projections.

Author: Casey W. Stark <caseywstark@gmail.com>
Copyright 2011.
Included September 6, 2011.

"""

import csv
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import yt
from yt.mods import *

from zeldovich import get_redshift, get_scale_factor, perturbed_x, perturbed_v_x

# path to initial conditions files
ic_path = os.path.abspath("../ics/")
run_path = os.path.abspath("../run_hydro/")

### Read in params
reader = csv.DictReader(open(os.path.join(ic_path, "zeldovich_params.csv"), "r"))
params = reader.__next__()

initial_z = float(params["initial_z"])
caustic_z = float(params["caustic_z"])
initial_a = get_scale_factor(initial_z)
box_length = float(params["box_length_long"])
box_length_short = float(params["box_length_short"])
k = float(params["k"])
offset = float(params["offset"])
H_0 = float(params["H_0"])
sheet_normal = np.array((params["normal_x"], params["normal_y"],
                                    params["normal_z"])).astype(float)

normal_length = min(box_length / sheet_normal[0],
                    box_length_short / sheet_normal[1],
                    box_length_short / sheet_normal[2])
q_array = np.linspace(0, normal_length, 100 * box_length * k)

print("")
print("Analyzing Zeldovich outputs.")

### Figure out plotfiles and figures
# If the plotfiles are in a different path, but the base directory here
pfs_path = run_path
fig_dir = os.path.abspath(os.path.dirname(__file__))

# Make a figs dir in the same dir as the pf's.
if not os.path.exists(os.path.join(fig_dir, "figs_hydro")):
    os.makedirs(os.path.join(fig_dir, "figs_hydro"))
figs_path = os.path.join(fig_dir, "figs_hydro")

# search for pltnnnnn dirs
pf_paths = [os.path.join(pfs_path, f) for f in os.listdir(pfs_path)
            if os.path.isdir(os.path.join(pfs_path, f)) and f.startswith("plt")]

# change dir to the figs dir so we can save there.
os.chdir(figs_path)

### Go through each particles file and test
for pf_path in pf_paths:
    plt_num = int(pf_path[-5:])

    print("")
    print("Processing pf{:05d}".format(plt_num))
    print("==================")
    print("")
    sys.stdout.write("Reading in particle data ... ")
    sys.stdout.flush()

    # read comoving_a
    a_file = open(os.path.join(pf_path, "comoving_a"), "r")
    current_a = float(a_file.readline())
    a_file.close()
    current_z = get_redshift(current_a)

    # open particles file
    p_file_path = os.path.join(pf_path, "particles.ascii")
    p_data = np.loadtxt(p_file_path, skiprows=5)

    print("DONE")
    sys.stdout.write("Processing particle data ... ")
    sys.stdout.flush()

    # reference columns into aliases and rotate.
    # note we divide the velocities by the scale factor to make them comoving.
    xs = p_data[:, 0]
    ys = p_data[:, 1]
    zs = p_data[:, 2]

    num_particles = len(xs)
    # Distances from the plane defined by the normal at (0, 0, 0)
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


    # get Zeldovich predictions
    zel_xs = perturbed_x(q_array + offset, current_z, caustic_z, k) - offset
    zel_vxs = perturbed_v_x(q_array + offset, current_z, caustic_z, k, H_0)

    print("DONE")
    print("")
    print("Checking velocities along the sheet ... ")

    max_v_along = np.max(v_alongs)

    print("Max velocity is {:.2e} km/s.".format(max_v_along))
    print("")
    sys.stdout.write("Making plots ... ")
    sys.stdout.flush()

    # Velocity histogram
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel("Velocity along sheet (km / s)")
    ax.set_ylabel("Number of particles")
    ax.hist(v_alongs)

    fig.savefig("velocity_along_{:05d}.png".format(plt_num))
    plt.close()

    # Phase space normal to sheets
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlim(0, box_length)
    ax.set_ylim(-1.2*zel_vxs.max(), 1.2*zel_vxs.max())
    ax.set_title("Zeldovich phase space")
    ax.set_xlabel("Comoving normal position (Mpc)")
    ax.set_ylabel("Comoving normal velocity (km/s)")
    ax.plot(pos_normals, v_normals, color="None", linestyle="None", marker=".",
            markerfacecolor=None, markeredgecolor="red", alpha=0.2)
    ax.plot(zel_xs, zel_vxs, color="gray", linestyle="-")
    ax.text(0.83, 0.92, r"$z_c =$ {:.2f}".format(caustic_z), transform=ax.transAxes, fontsize=16)
    ax.text(0.83, 0.85, r"$z =$ {:.2f}".format(current_z), transform=ax.transAxes, fontsize=16)

    fig.savefig("phase_space_{:05d}.png".format(plt_num))
    plt.close()

    # 3d particle positions
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_title("Zeldovich Particle Positions")
    ax.set_xlabel("x (Mpc)")
    ax.set_ylabel("y (Mpc)")
    ax.set_zlabel("z (Mpc)")
    ax.view_init(20, -85)
    ax.plot(xs, ys, zs, '.', alpha=0.015)
    ax.text(box_length, 0.0, box_length_short*1.25, r"$z_c =$ {:.2f}".format(caustic_z))
    ax.text(box_length, 0.0, box_length_short*1.1, r"$z =$ {:.2f}".format(current_z))
    ax.set_xlim3d(0, box_length)
    ax.set_ylim3d(0, box_length_short)
    ax.set_zlim3d(0, box_length_short)

    fig.savefig("positions_{:05d}.png".format(plt_num))
    plt.close()
    
    # Gas velocity, density, and temperature
    pf = load(pf_path)
    dims = pf.domain_dimensions
    dens = pf.covering_grid(0,left_edge=[0,0,0],dims=dims,fields=["density"])['density'].reshape((dims[0],dims[1],dims[2])).value
    mean_dens = np.mean(dens)
    temp = pf.covering_grid(0,left_edge=[0,0,0],dims=dims,fields=["Temp"])['Temp'].reshape((dims[0],dims[1],dims[2])).value
    xmom = pf.covering_grid(0,left_edge=[0,0,0],dims=dims,fields=["xmom"])['xmom'].reshape((dims[0],dims[1],dims[2])).value
    # Prepare the plot
    fig,ax = plt.subplots(3,1,figsize=(5,12))
    ax[0].set_xlim(0,box_length)
    ax[0].set_ylim(-1.2*(xmom/dens).max(),1.2*(xmom/dens).max())
    ax[0].set_title("Baryonic Zeldovich Pancake")
    ax[0].set_ylabel("Velocity [km/s]")
    ax[1].set_xlim(0,box_length)
    ax[1].set_ylim(0.75*dens.min()/mean_dens,1.5*dens.max()/mean_dens)
    ax[1].set_yscale("log")
    ax[1].set_ylabel("Overdensity")
    ax[2].set_ylabel("Temperature")
    ax[2].set_xlabel("Distance [cMpc]")
    ax[2].set_xlim(0,box_length)
    ax[2].set_ylim(0.75*temp.min(),1.5*temp.max())
    ax[2].set_yscale("log")
    # Loop over y/z columns?
    # Unnecessary because all columns are identical!
    #for ii in range(64):
    #    for jj in range(64):
    ax[0].plot((box_length/dims[0])*np.arange(dims[0]),np.roll(xmom[:,dims[1]//2-1,dims[2]//2-1]/dens[:,dims[1]//2-1,dims[2]//2-1],0), color="None", linestyle="None", marker=".",
           markerfacecolor=None, markeredgecolor="red", alpha=1.0)
    ax[1].plot((box_length/dims[0])*np.arange(dims[0]),np.roll(dens[:,dims[1]//2-1,dims[2]//2-1],0)/mean_dens, color="None", linestyle="None", marker=".",
           markerfacecolor=None, markeredgecolor="red", alpha=1.0,zorder=2)
    ax[2].plot((box_length/dims[0])*np.arange(dims[0]),np.roll(temp[:,dims[1]//2-1,dims[2]//2-1],0), color="None", linestyle="None", marker=".",
           markerfacecolor=None, markeredgecolor="red", alpha=1.0)

    ax[2].plot((box_length/dims[0])*np.arange(dims[0]),(184.1*((1+current_z)/101.0)**2)*(np.mean(dens,axis=(1,2))/mean_dens)**(2.0/3.0),
                   c='k',lw=1.5,zorder=1)
    fig.savefig("gas_{:05d}.png".format(plt_num))
    plt.close()
    
    print("DONE")

print("FULLY COMPLETE.")
