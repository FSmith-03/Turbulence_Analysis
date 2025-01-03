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


def main():
    # Constants
    L_e, tol, x_boundary, y_boundary, z_boundary, _, _, Nx, _, plot_limit, N_E = inp.constants()
    r = np.linspace(0, 2*x_boundary, 2*int(x_boundary/tol))
    r = r[:plot_limit]
    f_theoretical = np.exp(-r**2/L_e**2)
    g_theoretical = (1 - r**2/L_e**2) * np.exp(-r**2/L_e**2)
    f = f_theoretical
    g = g_theoretical
    return