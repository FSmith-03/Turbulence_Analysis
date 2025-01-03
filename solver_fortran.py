import input_settings as inp
import random_generator as rg
import positions_and_velocities_copy as pv
import numpy as np
import correlation_functions as cf
import matplotlib.pyplot as plt
import plot_correlations as pc
import energy_sepctrum as es
import fmodpy
import pandas as pd
import time

write = False

def write_to_csv(data, filename="output.csv"):
    # Convert each array in `data` to a DataFrame column
    df = pd.DataFrame(data)
    df.to_csv(filename, index=False)
    print(f"Data written to {filename}")

def filenamer(choice, test=None):
    x_boundary = inp.constants()[2]
    ARN = inp.constants()[-1]
    if test is None:
        if choice == 1:
            filename = f"velocitycomponents_{x_boundary}.csv"
        elif choice == 2:
            filename = f"correlationfunctions_{x_boundary}.csv"
    else:
        if choice == 1:
            filename = f"velocitycomponents_{x_boundary}_ARN_{ARN[test]}.csv"
        elif choice == 2:
            filename = f"correlationfunctions_{x_boundary}_ARN_{ARN[test]}.csv"
    return filename

def solver_needle():
    # Constants
    L_e, tol, x_boundary, y_boundary, z_boundary, _, _, Nx, _, plot_limit, N_E_base, ARN_list = inp.constants()
    for test in range(len(ARN_list)):
        print(f"Test number: {test}, ARN: {ARN_list[test]}")
        ARN = ARN_list[test]
        N_E = int(N_E_base / ARN)
        u_total = np.zeros(N_E)
        v_total = np.zeros(N_E)
        w_total = np.zeros(N_E)
        # Random angles and positions
        print("Total number of eddies: ", N_E)
        print("Generating random angles and positions")
        theta_list = rg.random_angles(N_E)
        a_list = rg.random_positions(x_boundary, y_boundary, z_boundary, N_E)
        print(a_list.shape)
        # Total velocities
        print("Calculating total velocities")
        start_time = time.time()  # Start the timer
        u_total, v_total, w_total = pv.total_velocities_2(Nx, x_boundary, N_E, ARN, test, tol, theta_list, a_list, "output.dat")
        elapsed_time = time.time() - start_time  # End the timer
        print(f"Time taken to calculate total velocities: {elapsed_time:.2f} seconds")
        # Write to csv
        if write == True:
            print("Writing data to csv")
            write_to_csv({"u_total": u_total, "v_total": v_total, "w_total": w_total,}, filenamer(1,test))
        #Mean velocity squared
        u_2_average = pv.u_2_average(u_total)
        v_2_average = pv.u_2_average(v_total)
        w_2_average = pv.u_2_average(w_total)
        print(u_2_average, v_2_average, w_2_average)
        # Correlation functions
        print("Calculating correlation functions")
        start_time = time.time()  # Start the timer
        r, f, g, f_s, max_index_plot = cf.correlation_functions_vect2(u_total, w_total, tol, plot_limit, x_boundary, u_2_average, w_2_average)
        elapsed_time = time.time() - start_time  # End the timer
        print(f"Time taken to calculate correlation functions: {elapsed_time:.2f} seconds")
        print("Calculating energy spectrum")
        f_filter, g_filter = es.f_and_g_filter(r, f, g)
        u_2_overall = (u_2_average + w_2_average) / 2
        E_k, k = es.energy_spectrum(r, f_filter, g_filter, tol, u_2_overall)
        total_energy = np.trapz(E_k, k)
        print(f"Integrated energy spectrum: {total_energy}")
        print(f"3/2 u^2: {3/2 * u_2_average**2}")
        print(f"Fractional difference: {((total_energy - 3/2 * u_2_average**2) / (3/2 * u_2_average**2))*100} %")
        es.energy_spectrum_plot(E_k, k, ARN)
        #Plot correlation and structure functions
        #Calculate Townsend and Signature functions
        print("Calculating Townsend and Signature functions")
        townsends, dvdr = cf.townsend_structure(f, r, u_2_average)
        signature_1 = cf.signature_function_1(f, r, u_2_average)
        print("Calculating Signature function 2")
        signature_2 = cf.signature_function_2(f, r, u_2_average, max_index_plot, x_boundary, tol)
        if write == True:
            print("Writing data to csv")
            E_k_write = np.zeros(len(r))
            k_write = np.zeros(len(r))
            E_k_write[:len(E_k)] = E_k
            k_write[:len(k)] = k
            write_to_csv({"r": r, "f": f, "g": g, "f_s": f_s, "E_k": E_k_write, "k": k_write, "V_t": townsends, "V_1": signature_1, "V_2": signature_2}, filenamer(2, test))
        print("Plotting correlation and structure functions")
        if test == 0:
            pc.theoretical_f(r, f, max_index_plot, L_e)
            pc.theoretical_g(r, g, max_index_plot, L_e)
        #pc.structure_plotter(r, f_s, max_index_plot, f'Structure Function (ARN={ARN_list[test]:.2f})')
        #pc.structure_plotter(r, townsends, max_index_plot, f'Townsend Structure Function (ARN={ARN_list[test]:.2f})')
        #pc.structure_plotter(r, signature_1, max_index_plot, f'Signature Function 1 (ARN={ARN_list[test]:.2f})')
        #pc.structure_plotter(r, signature_2, max_index_plot, f'Signature Function 2 (ARN={ARN_list[test]:.2f})')
    pc.plot_structure_comparison(max_index_plot, ARN_list, x_boundary, 'Townsend Structure Function')
    pc.plot_structure_comparison(max_index_plot, ARN_list, x_boundary, 'Signature Function 1')
    pc.plot_structure_comparison(max_index_plot, ARN_list, x_boundary, 'Signature Function 2')
    es.energy_spectrum_plot_comparison(ARN_list, max_index_plot, x_boundary)
    plt.show()
    return

def solver_sphere():
    # Constants
    L_e, tol, x_boundary, y_boundary, z_boundary, _, _, Nx, _, plot_limit, N_E_list, ARN_list = inp.constants()
    N_E = int(N_E_list)
    u_total = np.zeros(N_E)
    v_total = np.zeros(N_E)
    w_total = np.zeros(N_E)
    # Random angles and positions
    print("Total number of eddies: ", N_E)
    print("Generating random angles and positions")
    theta_list = rg.random_angles(N_E)
    a_list = rg.random_positions(x_boundary, y_boundary, z_boundary, N_E)
    print(a_list.shape)
    # Total velocities
    print("Calculating total velocities")
    start_time = time.time()  # Start the timer
    u_total, v_total, w_total = pv.total_velocities(Nx, x_boundary, tol, N_E, theta_list, a_list)
    elapsed_time = time.time() - start_time  # End the timer
    print(f"Time taken to calculate total velocities: {elapsed_time:.2f} seconds")
    # Write to csv
    print("Writing data to csv")
    write_to_csv({"u_total": u_total, "v_total": v_total, "w_total": w_total,}, filenamer(1))
    #Mean velocity squared
    u_2_average = pv.u_2_average(u_total)
    v_2_average = pv.u_2_average(v_total)
    w_2_average = pv.u_2_average(w_total)
    print(u_2_average, v_2_average, w_2_average)
    # Correlation functions
    print("Calculating correlation functions")
    start_time = time.time()  # Start the timer
    r, f, g, f_s, max_index_plot = cf.correlation_functions_vect2(u_total, v_total, tol, plot_limit, x_boundary, u_2_average, v_2_average)
    elapsed_time = time.time() - start_time  # End the timer
    print(f"Time taken to calculate correlation functions: {elapsed_time:.2f} seconds")
    print("Calculating energy spectrum")
    f_filter, g_filter = es.f_and_g_filter(r, f, g)
    E_k, k = es.energy_spectrum(r, f_filter, g_filter, tol, u_2_average)
    total_energy = np.trapz(E_k, k)
    print(f"Integrated energy spectrum: {total_energy}")
    print(f"3/2 u^2: {3/2 * u_2_average**2}")
    print(f"Fractional difference: {((total_energy - 3/2 * u_2_average**2) / (3/2 * u_2_average**2))*100} %")
    es.energy_spectrum_plot(E_k, k)
    #Plot correlation and structure functions
    #Calculate Townsend and Signature functions
    print("Calculating Townsend and Signature functions")
    townsends, dvdr = cf.townsend_structure(f, r, u_2_average)
    signature_1 = cf.signature_function_1(f, r, u_2_average)
    print("Calculating Signature function 2")
    signature_2 = cf.signature_function_2(f, r, u_2_average, max_index_plot, x_boundary, tol)
    print("Writing data to csv")
    E_k_write = np.zeros(len(r))
    k_write = np.zeros(len(r))
    E_k_write[:len(E_k)] = E_k
    k_write[:len(k)] = k
    write_to_csv({"r": r, "f": f, "g": g, "f_s": f_s, "E_k": E_k_write, "k": k_write, "V_t": townsends, "V_1": signature_1, "V_2": signature_2}, filenamer(2))
    print("Plotting correlation and structure functions")
    pc.theoretical_f(r, f, max_index_plot, L_e)
    pc.theoretical_g(r, g, max_index_plot, L_e)
    pc.structure_plotter(r, f_s, max_index_plot, 'Structure Function')
    pc.structure_plotter(r, townsends, max_index_plot, 'Townsend Structure Function')
    pc.structure_plotter(r, signature_1, max_index_plot, 'Signature Function 1')
    pc.structure_plotter(r, signature_2, max_index_plot, 'Signature Function 2')
    plt.show()
    return

def processer():
    L_e, tol, x_boundary, y_boundary, z_boundary, _, _, Nx, _, plot_limit, N_E_list, ARN_list = inp.constants()
    velocities_total = np.loadtxt("output.dat")
    print(len(velocities_total))
    u_total = velocities_total[:Nx]
    v_total = velocities_total[Nx:2*Nx]
    w_total = velocities_total[2*Nx:]
    print("Writing data to csv")
    write_to_csv({"u_total": u_total, "v_total": v_total, "w_total": w_total,}, filenamer(1))
    #Mean velocity squared
    u_2_average = pv.u_2_average(u_total)
    v_2_average = pv.u_2_average(v_total)
    w_2_average = pv.u_2_average(w_total)
    print(u_2_average, v_2_average, w_2_average)
    # Correlation functions
    print("Calculating correlation functions")
    start_time = time.time()  # Start the timer
    r, f, g, f_s, max_index_plot = cf.correlation_functions_vect2(u_total, w_total, tol, plot_limit, x_boundary, u_2_average, w_2_average)
    elapsed_time = time.time() - start_time  # End the timer
    print(f"Time taken to calculate correlation functions: {elapsed_time:.2f} seconds")
    # Write to csv
    print("Writing data to csv")
    write_to_csv({"r": r, "f": f, "g": g, "f_s": f_s}, filenamer(2))
    print("Calculating energy spectrum")
    f_filter, g_filter = es.f_and_g_filter(r, f, g)
    E_k, k = es.energy_spectrum(r, f_filter, g_filter, tol, u_2_average)
    total_energy = np.trapz(E_k, k)
    print(f"Integrated energy spectrum: {total_energy}")
    print(f"3/2 u^2: {3/2 * u_2_average**2}")
    print(f"Fractional difference: {((total_energy - 3/2 * u_2_average**2) / (3/2 * u_2_average**2))*100} %")
    es.energy_spectrum_plot(E_k, k)
    #Plot correlation and structure functions
    #Calculate Townsend and Signature functions
    print("Calculating Townsend and Signature functions")
    townsends, dvdr = cf.townsend_structure(f, r, u_2_average)
    signature_1 = cf.signature_function_1(f, r, u_2_average)
    print("Calculating Signature function 2")
    signature_2 = cf.signature_function_2(f, r, u_2_average, max_index_plot, x_boundary, tol)
    print("Plotting correlation and structure functions")
    pc.theoretical_f(r, f, max_index_plot, L_e)
    pc.theoretical_g(r, g, max_index_plot, L_e)
    pc.structure_plotter(r, f_s, max_index_plot, 'Structure Function')
    pc.structure_plotter(r, townsends, max_index_plot, 'Townsend Structure Function')
    pc.structure_plotter(r, signature_1, max_index_plot, 'Signature Function 1')
    pc.structure_plotter(r, signature_2, max_index_plot, 'Signature Function 2')
    plt.show()
    return

solver_needle()
#processer()