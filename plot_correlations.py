import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def structure_plotter(r, f, max_index, title):
    f = f[0:max_index]
    r = r[0:max_index]
    max_location = r[np.argmax(f)]
    max_location = round(max_location, 4)
    fig, ax = plt.subplots()
    ax.plot(r, f)
    ax.axvline(x=max_location, color='red', linestyle='--', label=f'Location of max value: {max_location}')
    ax.set_xlabel('r/L')
    ax.set_ylabel(title)
    plt.title(title)
    plt.legend()

def theoretical_f(r, f, max_index, L_e):
    r = r[0:max_index]
    f = f[0:max_index]
    f_theoretical = np.exp(-r**2/L_e**2)
    fig, ax = plt.subplots()
    ax.plot(r, f, label='f', c='r')
    ax.plot(r, f_theoretical, label='Theoretical f', c = 'b', linestyle='--')
    ax.set_xlabel('r/L_e')
    ax.set_ylabel('f')
    plt.legend()
    plt.title('Longitudinal Correlation Function')

def theoretical_g(r, g, max_index, L_e):
    g = g[0:max_index]
    r = r[0:max_index]
    g_theoretical = (1 - r**2/L_e**2) * np.exp(-r**2/L_e**2)
    fig, ax = plt.subplots()
    ax.plot(r, g, label='g', c='r')
    ax.plot(r, g_theoretical, label='Theoretical g', c = 'b', linestyle='--')
    ax.set_xlabel('r/L_e')
    ax.set_ylabel('g')
    plt.legend()
    plt.title('Transverse Correlation Function')

def plot_structure_comparison(max_index, ARN_list, x_boundary, title):
    fig, ax = plt.subplots()
    allocated = False
    
    for ARN in ARN_list:
        file = f"correlationfunctions_{x_boundary}_ARN_{ARN}.csv"
        data = pd.read_csv(file)

        if allocated == False:
            r = np.array(data["r"].values)
            r = r[0:max_index]
            allocated = True

        if title == "Townsend Structure Function":
            V_t = np.array(data["V_t"].values)
            V_t = V_t[0:max_index]
            ax.plot(r, V_t, label=f"ARN = {ARN}")
        
        if title == "Signature Function 1":
            V_1 = np.array(data["V_1"].values)
            V_1 = V_1[0:max_index]
            ax.plot(r, V_1, label=f"ARN = {ARN}")
        
        if title == "Signature Function 2":
            V_2 = np.array(data["V_2"].values)
            V_2 = V_2[0:max_index]
            ax.plot(r, V_2, label=f"ARN = {ARN}")
        
    ax.set_xlabel('r/L_e')
    ax.set_ylabel(title)
    plt.legend()
    plt.title(title)

        
    
    


            
        