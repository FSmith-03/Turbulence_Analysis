import numpy as np
import correlation_functions as cf
import matplotlib.pyplot as plt

def main():
    L_e = 1.0
    r = np.linspace(0, 5*L_e, int(5/0.05))
    f_theoretical = np.exp(-r**2/L_e**2)
    g_theoretical = (1 - r**2/L_e**2) * np.exp(-r**2/L_e**2)
    f_real = f_theoretical + np.random.normal(0, 0.0002*max(f_theoretical), len(r))
    g_real = g_theoretical + np.random.normal(0, 0.2, len(r))
    signature = cf.signature_function_1(f_theoretical, r, 1)
    signature_real = cf.signature_function_1(f_real, r, 1)
    fig, ax = plt.subplots()
    ax.plot(r, signature, label='Signature Function with no noise')
    ax.plot(r, signature_real, label='Signature Function with noise')
    ax.set_xlabel('r/L_e')
    ax.set_ylabel('Signature Function')
    plt.legend()
    plt.title('Comparison of Signature Functions with and without noise')
    fig, ax = plt.subplots()
    ax.plot(r, f_theoretical, label='Theoretical f')
    ax.plot(r, f_real, label='Real f')
    ax.set_xlabel('r/L_e')
    ax.set_ylabel('f')
    plt.legend()
    plt.title('Comparison of Theoretical and Real f')
    plt.show()
main()