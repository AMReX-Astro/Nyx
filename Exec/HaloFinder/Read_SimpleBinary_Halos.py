import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def read_binary_points(filename):
	with open(filename, 'rb') as file:
		# Read the binary data directly, assuming big-endian float32
		data = np.fromfile(file, dtype='>f4')  # '>f4' indicates big-endian float32

		# Reshape the data to Nx3 (x, y, z for each point)
		if len(data) % 5 != 0:
			raise ValueError("Data size is not a multiple of 5. The file might be corrupted or improperly formatted.")

		points = data.reshape((-1, 5))

		num_points = len(data) // 5
		print(f"Total number of points: {num_points}")

	return points, num_points

def plot_points(x, y, z):
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection='3d')

    # Scatter plot for points
    ax.scatter(x, y, z, c='blue', marker='o', s=1, alpha=0.8)

    # Set labels for the axes
    ax.set_xlabel('X Coordinate')
    ax.set_ylabel('Y Coordinate')
    ax.set_zlabel('Z Coordinate')

    ax.set_title('Lightcone shell')
    plt.show()


# Example usage
filename = "./Output/Halos/SimpleBinary/reeber_halos_0000177.bin"
points, num_points = read_binary_points(filename)

# Display the first few points
x = np.zeros(num_points);
y = np.zeros(num_points);
z = np.zeros(num_points);
xcen = 3850.0;
ycen = 3850.0;
zcen = 3850.0;
for i, point in enumerate(points[:num_points]):  # Adjust the slice as needed
	x[i], y[i], z[i], mass, n_cells = point
	#rad = np.sqrt((x[i] - xcen) * (x[i] - xcen) + (y[i] - ycen) * (y[i] - ycen) + (z[i] - zcen) * (z[i] - zcen));
	print(f"{x[i]:.15g}, {y[i]:.15g}, {z[i]:.15g}, {mass:.15g}, {n_cells:.15g}")
	#print(f"{x[i]:.15g}, {y[i]:.15g}, {z[i]:.15g}, {vx:.15g}, {vy:.15g}, {vz:.15g}")

#plot_points(x, y, z)


