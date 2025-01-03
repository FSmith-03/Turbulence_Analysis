import input_settings as inp
import random_generator as rg
import positions_and_velocities as pv
import numpy as np
import correlation_functions as cf
import matplotlib.pyplot as plt
import plot_correlations as pc
def solver():
    # Constants
    L_e, tol, x_boundary, y_boundary, z_boundary, _, _, Nx, _, plot_limit, N_E = inp.constants()
    # Random angles and positions
    print("Total number of eddies: ", N_E)
    print("Generating random angles and positions")
    theta_list = rg.random_angles(N_E)
    a_list = rg.random_positions(x_boundary, y_boundary, z_boundary, N_E)
    # Total velocities
    print("Calculating total velocities")
    u_total, v_total, w_total = pv.total_velocities(Nx, x_boundary, N_E, theta_list, a_list)
    #Mean velocity squared
    u_2_average = pv.u_2_average(u_total)
    v_2_average = pv.u_2_average(v_total)
    w_2_average = pv.u_2_average(w_total)
    print(u_2_average, v_2_average, w_2_average)
    # Correlation functions
    print("Calculating correlation functions")
    r, f, g, f_s, max_index_plot = cf.correlation_functions_vect(u_total, v_total, tol, plot_limit, u_2_average, v_2_average)
    #Plot correlation and structure functions
    print("Plotting correlation and structure functions")
    pc.theoretical_f(r, f, max_index_plot, L_e)
    pc.theoretical_g(r, g, max_index_plot, L_e)
    pc.structure_plotter(r, f_s, max_index_plot)
    plt.show()
    return

#solver()

def checks():
    L_e, tol, x_boundary, y_boundary, z_boundary, _, _, Nx, _, plot_limit, N_E = inp.constants()
    theta_list = rg.random_angles(N_E)
    a_list = rg.random_positions(x_boundary, y_boundary, z_boundary, N_E)
    pv.rotation_check(Nx, x_boundary, N_E, theta_list)

checks()