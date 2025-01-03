import numpy as np
from numba import njit
import sympy as sp

# Function to calculate the correlation functions
# Note only need u^2 average as the turbulence is isotropic

def correlation_functions_vect(u_total, v_total, tol, plot_limit, x_boundary, u_2_average, v_2_average):
    N_u = len(u_total)
    s = np.arange(N_u)
    max_index = int(plot_limit/tol)
    r = np.linspace(0, 2*x_boundary, 2*int(x_boundary/tol))
    f = np.zeros(N_u)
    g = np.zeros(N_u)
    f_s = np.zeros(N_u)
    
    for i in range(N_u):
        shift = s[i]
        if shift < N_u:
            u_shifted = u_total[shift:]
            v_shifted = v_total[shift:]
            valid_range = slice(0, N_u - shift)
            product_list_f = u_total[valid_range] * u_shifted[valid_range]
            product_list_g = v_total[valid_range] * v_shifted[valid_range]
            product_list_f_s = (u_total[valid_range] - u_shifted[valid_range])**2
            
            f[i] = np.mean(product_list_f)
            g[i] = np.mean(product_list_g)
            f_s[i] = np.mean(product_list_f_s)
    
    # Normalize the correlation functions
    f = f / u_2_average
    g = g / v_2_average
    return r, f, g, f_s, max_index

@njit
def correlation_functions_vect2(u_total, v_total, tol, plot_limit, x_boundary, u_2_average, v_2_average):
    N_u = len(u_total)
    max_index = int(plot_limit / tol)
    r = np.linspace(0, 2 * x_boundary, 2 * int(x_boundary / tol))
    
    f = np.zeros(N_u)
    g = np.zeros(N_u)
    f_s = np.zeros(N_u)
    
    
    for shift in range(N_u):
        # Calculate the valid range
        valid_length = N_u - shift
        
        # Accumulators for sums
        sum_f = 0.0
        sum_g = 0.0
        sum_f_s = 0.0
        
        # Perform element-wise operations manually for the valid range
        for i in range(valid_length):
            sum_f += u_total[i] * u_total[i + shift]
            sum_g += v_total[i] * v_total[i + shift]
            sum_f_s += (u_total[i] - u_total[i + shift])**2
        
        # Compute the means
        f[shift] = sum_f / valid_length
        g[shift] = sum_g / valid_length
        f_s[shift] = sum_f_s / valid_length
    
    # Normalize the correlation functions
    f /= u_2_average
    g /= v_2_average
    
    return r, f, g, f_s, max_index

def correlation_functions_base(u_total, v_total, tol, plot_limit, x_boundary, u_2_average, v_2_average):
    N_u = len(u_total)
    s = np.arange(N_u)
    max_index = int(plot_limit/tol)
    r = np.linspace(0, 2*x_boundary, 2*int(x_boundary/tol))
    f = []
    g = []
    f_s = []
    
    for i in range(N_u):
        shift = s[i]
        if shift < N_u:
            u_shifted = np.roll(u_total, -shift)
            v_shifted = np.roll(v_total, -shift)
            valid_range = slice(0, N_u - shift)
            product_list_f = u_total[valid_range] * u_shifted[valid_range]
            product_list_g = v_total[valid_range] * v_shifted[valid_range]
            product_list_f_s = (u_total[valid_range] - u_shifted[valid_range])**2
            
            f.append(np.mean(product_list_f))
            g.append(np.mean(product_list_g))
            f_s.append(np.mean(product_list_f_s))
    
    # Normalize the correlation functions
    f = f / u_2_average
    g = g / v_2_average
    return r, f, g, f_s, max_index


# Townsend's structure function
def townsend_structure(f, r, u_2_average):
    delta_v = 2 * u_2_average * (1 - f)
    dvdr = np.gradient(delta_v, r)
    townsends = 3/4 * dvdr
    return townsends, dvdr

def signature_function_1(f, r, u_2_average):
    signature_1 = np.zeros(len(r))
    dfdr = np.gradient(f, r)
    df2dr2 = np.gradient(dfdr, r)
    signature_1[1:] = 3/4 * u_2_average * (r[1:] * df2dr2[1:] - dfdr[1:])
    return signature_1

@njit
def g(x):
    return (7 * x**5 - 2 * x**7) * np.exp(-x**2)

@njit
def dgds(s, r):
    return ((35 / r**5) * s**4 - (14 / r**7) * s**6) * np.exp(-(s / r)**2) - (2 * s / r**2) * g(s / r)

@njit
def signature_function_2(f, r, u_2_average, limit, x_boundary, tol):
    # Use a slice of r directly, avoid reassigning
    r = np.linspace(0, 2*x_boundary, 2*int(x_boundary/tol))
    s = r
    n_r = len(s)
    n_s = len(s)

    delta_v_2 = 2 * u_2_average * (1 - f)

    # Initialize arrays
    u_v = np.zeros(n_r)
    integrand = np.zeros(n_r)
    signature_2 = np.zeros(n_r)

    # Compute u_v and integrand
    for j in range(limit):
        r_j = s[j]
        if r_j == 0:
            continue  # Skip if r_j is 0

        # Compute u_v for the current r_j
        u_v[j] = delta_v_2[-1] * g(s[-1] / r_j) - delta_v_2[0] * g(s[0] / r_j)

        # Compute the integrand for all s
        for i in range(n_s):
            integrand[i] = dgds(s[i], r_j) * delta_v_2[i]
        
        # Compute the signature function
        signature_2[j+1] = u_v[j] - np.trapz(integrand[1:], s[1:])

    return signature_2

@njit
def integral_by_parts(u, du_dx, v, x, max_plot_index):
    """
    Perform numerical integration by parts with manual trapezoidal rule (Compatible with njit).

    Parameters:
    - u: array-like, u(x) function values
    - du_dx: array-like, derivative of u(x) with respect to x
    - v: array-like, v(x) function values
    - x: array-like, x values for numerical integration

    Returns:
    - result: float, value of the integral
    """
    # Boundary term: [u * v]
    boundary_term = u[-1] * v[-1] - u[0] * v[0]  # Swap order due to limits

    # Trapezoidal rule for the integral of u' * v
    integral_term = 0.0
    for i in range(1, max_plot_index):
        dx = x[i] - x[i - 1]
        integral_term += (du_dx[i] * v[i] + du_dx[i - 1] * v[i - 1]) * 0.5 * dx

    # Integration by parts formula: [u * v] - ∫ u' * v dx
    return boundary_term - integral_term

@njit
def compute_gradient(f, s):
    """
    Compute the gradient df/ds using central differences for the interior points
    and forward/backward differences for the boundaries.
    
    Parameters:
        f (numpy.ndarray): Array of function values.
        s (numpy.ndarray): Array of independent variable values.
    
    Returns:
        numpy.ndarray: Gradient df/ds.
    """
    n = len(f)
    df_ds = np.zeros_like(f)
    for i in range(n):
        if i == 0:
            # Forward difference at the first point
            df_ds[i] = (f[i + 1] - f[i]) / (s[i + 1] - s[i])
        elif i == n - 1:
            # Backward difference at the last point
            df_ds[i] = (f[i] - f[i - 1]) / (s[i] - s[i - 1])
        else:
            # Central difference for interior points
            df_ds[i] = (f[i + 1] - f[i - 1]) / (s[i + 1] - s[i - 1])
    return df_ds

@njit
def compute_integral(f, s):
    """
    Compute the integral of f with respect to s using the trapezoidal rule.
    
    Parameters:
        f (numpy.ndarray): Array of function values.
        s (numpy.ndarray): Array of independent variable values.
    
    Returns:
        float: Integral of f with respect to s.
    """
    n = len(f)
    integral = 0.0
    for i in range(1, n):
        integral += 0.5 * (f[i] + f[i - 1]) * (s[i] - s[i - 1])
    return integral

@njit
def signature_function_2_1(delta_v_2, r, u_2_average, limit):
    # Use a slice of r directly, avoid reassigning
    r = r[1:]
    s = r
    n_r = len(s)
    signature_2 = np.zeros(n_r+1)
    for i,r_i in enumerate(r):
        g_v = g(s/r_i)
        dvds = compute_gradient(delta_v_2, s)
        integrand = g_v * dvds
        signature_2[i+1] = compute_integral(integrand, s) / r_i
    return signature_2
