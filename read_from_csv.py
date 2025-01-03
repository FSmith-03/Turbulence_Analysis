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

def read_and_plot(N):
    L_e = inp.constants()[0]
    filename_1 = f"velocitycomponents_{N}.csv"
    filename_2 = f"correlationfunctions_{N}.csv"
    data_1 = pd.read_csv(filename_1)
    data_2 = pd.read_csv(filename_2)
    u_total = data_1["u_total"]
    v_total = data_1["v_total"]
    w_total = data_1["w_total"]
    r = data_2["r"]
    f = data_2["f"]
    g = data_2["g"]
    f_s = data_2["f_s"]
    # Convert all data to numpy arrays
    u_total = np.array(u_total)
    v_total = np.array(v_total)
    w_total = np.array(w_total)
    r = np.array(r)
    f = np.array(f)
    g = np.array(g)
    f_s = np.array(f_s)
    max_index_plot = len(r)
    #Mean velocity squared
    u_2_average = pv.u_2_average(u_total)
    v_2_average = pv.u_2_average(v_total)
    w_2_average = pv.u_2_average(w_total)
    print(u_2_average, v_2_average, w_2_average)
    # Filter the correlation functions
    f_filter, g_filter = es.f_and_g_filter(r, f, g)
    f = f_filter
    g = g_filter
    #Calculate Townsend and Signature functions
    print("Calculating Townsend and Signature functions")
    townsends, dvdr = cf.townsend_structure(f, r, u_2_average)
    signature_1 = cf.signature_function_1(f, r, u_2_average)
    signature_2 = cf.signature_function_2(f_s, r, u_2_average, max_index_plot)
    print("Plotting correlation and structure functions")
    pc.theoretical_f(r, f, max_index_plot, L_e)
    pc.theoretical_g(r, g, max_index_plot, L_e)
    pc.structure_plotter(r, f_s, max_index_plot, 'Structure Function')
    pc.structure_plotter(r, townsends, max_index_plot, 'Townsend Structure Function')
    pc.structure_plotter(r, signature_1, max_index_plot, 'Signature Function 1')
    pc.structure_plotter(r, signature_2, max_index_plot, 'Signature Function 2')
    plt.show()
    return

read_and_plot(1100)