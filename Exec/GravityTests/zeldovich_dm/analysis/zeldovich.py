import numpy as np

def perturbed_x(q, z, z_c, k):
    return (q - (1.0 + z_c) / (1.0 + z) * np.sin(k * q) / k)

def perturbed_v_x(q, z, z_c, k, H_0):
    return (-H_0 * (1.0 + z_c) * np.sqrt(1.0 + z) * np.sin(k * q) / (k))

def get_scale_factor(z):
    return 1.0 / (1.0 + z)

def get_redshift(a):
    return 1.0 / a - 1.0
