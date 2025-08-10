import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys
import os

def read_binary_points(filename):
    with open(filename, 'rb') as file:
        data = np.fromfile(file, dtype='>f4')  # big-endian float32

        if len(data) % 5 != 0:
            raise ValueError("Data size is not a multiple of 5. The file might be corrupted or improperly formatted.")

        points = data.reshape((-1, 5))
        num_points = len(points)
        print(f"Total number of points: {num_points}")

    return points, num_points

def plot_points(x, y, z):
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x, y, z, c='blue', marker='o', s=1, alpha=0.8)
    ax.set_xlabel('X Coordinate')
    ax.set_ylabel('Y Coordinate')
    ax.set_zlabel('Z Coordinate')
    ax.set_title('Lightcone Shell')
    plt.show()

def main():
    if len(sys.argv) != 2:
        print(f"Usage: python {os.path.basename(sys.argv[0])} <binary_file>")
        sys.exit(1)

    filename = sys.argv[1]

    if not os.path.isfile(filename):
        print(f"Error: File '{filename}' not found.")
        sys.exit(1)

    points, num_points = read_binary_points(filename)

    # Extract coordinates
    x = np.zeros(num_points)
    y = np.zeros(num_points)
    z = np.zeros(num_points)

    for i, point in enumerate(points):
        x[i], y[i], z[i], mass, n_cells = point
        print(f"{x[i]:.15g}, {y[i]:.15g}, {z[i]:.15g}, {mass:.15g}, {n_cells:.15g}")

    plot_points(x, y, z)

if __name__ == "__main__":
    main()
