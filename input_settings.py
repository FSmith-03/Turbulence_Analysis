# This file contains the input settings for the simulation. 
# The user can change the values of the variables to suit their needs. 
# The variables are then used in the main file to run the simulation.

# Function to calculate the number of eddies to fill the box with one
# eddy per unit volume
import numpy as np

def spacefiller(x_boundary, y_boundary, z_boundary, ARN=None):
    N = 8*x_boundary*y_boundary*z_boundary
    return N

def constants():
    ARN = [1, 1.5, 2, 5, 10, 20, 40]
    # ARN = [1]
    L_e = 1
    tol = 0.05
    #Change x length to different values to change the size of the box
    # y and z set to 2.95
    x_boundary = 500
    y_boundary = 2.95
    z_boundary = 2.95
    Nxf = 2*x_boundary/tol
    Nyz = 2*y_boundary/tol
    #Convert Nxf to an integer to be used in an array input
    Nx = int(Nxf)
    Nyz = int(Nyz)
    limit = 5
    N_E = spacefiller(x_boundary, y_boundary, z_boundary)
    return L_e, tol, x_boundary, y_boundary, z_boundary, Nxf, Nyz, Nx, Nyz, limit, N_E, ARN