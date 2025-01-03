import numpy as np
import rotation_matrix as rm
import numba as nb
import random_generator as rg
from joblib import Parallel, delayed
import fmodpy
import matplotlib.pyplot as plt
import ctypes

#print(dir(code))
def sensor_line(x_boundary, Nxf):
    x = np.linspace(0, x_boundary, Nxf)
    y = np.zeros_like(x)
    z = np.zeros_like(x)
    pos_vector = np.vstack([x, y, z])
    return pos_vector

def eddy_range(pos_vectors):
    factor = pos_vectors[0]**2 + pos_vectors[1]**2 + pos_vectors[2]**2
    mask = np.sqrt(factor) < 4
    first_index = np.argmax(mask)
    last_index = len(mask) - np.argmax(mask[::-1])
    return first_index, last_index, factor[first_index:last_index]

def velocity_generator(xaxis_trimmed, factor):
    u = -xaxis_trimmed[1] * np.exp(-factor * 2)
    v = xaxis_trimmed[0] * np.exp(-factor * 2)
    w = np.zeros_like(u)
    return u, v, w

def u_2_average(u_total):
    u_2 = u_total**2
    u_2_average = np.mean(u_2)
    return u_2_average

def total_velocities(Nx, x_boundary, tol, N_E, theta_list, a_list):
    code = fmodpy.fimport("routines_sphere.f90")
    velocities_total = np.zeros((Nx, 3))
    xaxis = sensor_line(x_boundary, Nx)
    L_e = 1
    ARP = 1
    input = np.array([x_boundary, Nx, N_E])
    input = np.asarray(input, dtype=ctypes.c_int, order='F')
    a_list = np.asarray(a_list, dtype=ctypes.c_double, order='F')
    theta_list = np.asarray(theta_list, dtype=ctypes.c_double, order='F')
    code.main_calculation(input, a_list, theta_list)
    velocities_total = np.loadtxt("output.dat")
    print("Shape of velocities_total:", np.shape(velocities_total))

    # Separate components
    u_total = velocities_total[:Nx]
    v_total = velocities_total[Nx:2*Nx]
    w_total = velocities_total[2*Nx:]
    #u_total = velocities_total[:,0]
    #v_total = velocities_total[:,1]
    #w_total = velocities_total[:,2]
    #plot u v and w
    plt.plot(xaxis[0], u_total, label='u')
    #plt.plot(xaxis[0], v_total, label='v')
    #plt.plot(xaxis[0], w_total, label='w')
    plt.legend()
    return u_total, v_total, w_total

def total_velocities_2(Nx, x_boundary, N_E, ARN, test, tol, theta_list, a_list, filename):
    code = fmodpy.fimport("routines_needles.f90")
    xaxis = sensor_line(x_boundary, Nx)
    filename = "output.dat"
    L_e = 1
    ARP = 1
    input = np.array([x_boundary, Nx, int(N_E), int(ARP), int(ARN)])
    input = np.asarray(input, dtype=ctypes.c_int, order='F')
    a_list = np.asarray(a_list, dtype=ctypes.c_float, order='F')
    theta_list = np.asarray(theta_list, dtype=ctypes.c_float, order='F')
    code.main_calculation(input, a_list, theta_list)
    print("filename:", filename)
    velocities_total = np.loadtxt(filename)
    print("Shape of velocities_total:", np.shape(velocities_total))

    # Separate components
    print(Nx)
    u_total = velocities_total[:Nx]
    v_total = velocities_total[Nx:2*Nx]
    w_total = velocities_total[2*Nx:]
    # plot u v and w
    #plt.plot(xaxis[0], u_total, label='u')
    #plt.plot(xaxis[0], v_total, label='v')
    #plt.plot(xaxis[0], w_total, label='w')
    #plt.legend()
    return u_total, v_total, w_total

def process_single_eddy(Nx, x_boundary, theta_list, a_list, i):
    u_result, v_result, w_result = np.zeros(Nx), np.zeros(Nx), np.zeros(Nx)
    xaxis = sensor_line(x_boundary, Nx)
    a  = a_list[:,i][:, np.newaxis]
    eddy_pos_translated = xaxis - a
    first_index, last_index, factor = eddy_range(eddy_pos_translated)
    eddy_pos_trimmed = eddy_pos_translated[:, first_index:last_index]
    thetax, thetay, thetaz = theta_list[i]
    R = rm.rotation_total(thetax, thetay, thetaz)
    R_inv = rm.rotation_total(-thetax, -thetay, -thetaz)
    eddy_pos_rotated = R @ eddy_pos_trimmed
    u, v, w = velocity_generator(eddy_pos_rotated, factor)
    u_rotated, v_rotated, w_rotated = R_inv @ np.array([u, v, w])
    u_result[first_index:last_index] = u_rotated
    v_result[first_index:last_index] = v_rotated
    w_result[first_index:last_index] = w_rotated
    return u_result, v_result, w_result

