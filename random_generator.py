import random
import numpy as np


def random_angles(N_eddies):
    theta_x = np.random.uniform(0, 2 * np.pi, (N_eddies))
    theta_y = np.random.uniform(0, 2 * np.pi, (N_eddies))
    theta_z = np.random.uniform(0, 2 * np.pi, (N_eddies))
    angles_list = np.array([theta_x, theta_y, theta_z])
    return angles_list


def random_positions(x_boundary, y_boundary, z_boundary, N_eddies):
    positions_list_x = np.random.uniform(-x_boundary, x_boundary, (N_eddies))
    positions_list_y = np.random.uniform(-y_boundary, y_boundary, (N_eddies))
    positions_list_z = np.random.uniform(-z_boundary, z_boundary, (N_eddies))
    positions_list = np.array([positions_list_x, positions_list_y, positions_list_z])
    return positions_list
