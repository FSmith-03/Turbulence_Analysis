import numpy as np
import pandas as pd
import positions_and_velocities as pv
import correlation_functions as cf
import plot_correlations as pc
import input_settings as inp
import matplotlib.pyplot as plt
def load_column_from_csv(file_name, column_name):
    data = pd.read_csv(file_name)
    entry = data[column_name].values
    entry = np.array(entry)
    return entry

u_total = load_column_from_csv("output1_20000.csv", "u_total")
v_total = load_column_from_csv("output1_20000.csv", "v_total")
w_total = load_column_from_csv("output1_20000.csv", "w_total")
u_2_average = pv.u_2_average(u_total)
v_2_average = pv.u_2_average(v_total)
w_2_average = pv.u_2_average(w_total)
L_e, tol, x_boundary, y_boundary, z_boundary, _, _, Nx, _, plot_limit, N_E = inp.constants()
r, f, g, f_s, max_index = cf.correlation_functions_vect(u_total, v_total, tol, plot_limit, u_2_average, v_2_average)
pc.theoretical_f(r, f, max_index, L_e)
pc.theoretical_g(r, g, max_index, L_e)
pc.structure_plotter(r, f_s, max_index)
plt.show()