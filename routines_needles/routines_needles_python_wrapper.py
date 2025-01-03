'''This Python code is an automatically generated wrapper
for Fortran code made by 'fmodpy'. The original documentation
for the Fortran source code follows.


'''

import os
import ctypes
import platform
import numpy

# --------------------------------------------------------------------
#               CONFIGURATION
# 
_verbose = True
_fort_compiler = "gfortran"
_shared_object_name = "routines_needles." + platform.machine() + ".so"
_this_directory = os.path.dirname(os.path.abspath(__file__))
_path_to_lib = os.path.join(_this_directory, _shared_object_name)
_compile_options = ['-fPIC', '-shared', '-O3']
_ordered_dependencies = ['routines_needles.f90', 'routines_needles_c_wrapper.f90']
_symbol_files = []# 
# --------------------------------------------------------------------
#               AUTO-COMPILING
#
# Try to import the prerequisite symbols for the compiled code.
for _ in _symbol_files:
    _ = ctypes.CDLL(os.path.join(_this_directory, _), mode=ctypes.RTLD_GLOBAL)
# Try to import the existing object. If that fails, recompile and then try.
try:
    # Check to see if the source files have been modified and a recompilation is needed.
    if (max(max([0]+[os.path.getmtime(os.path.realpath(os.path.join(_this_directory,_))) for _ in _symbol_files]),
            max([0]+[os.path.getmtime(os.path.realpath(os.path.join(_this_directory,_))) for _ in _ordered_dependencies]))
        > os.path.getmtime(_path_to_lib)):
        print()
        print("WARNING: Recompiling because the modification time of a source file is newer than the library.", flush=True)
        print()
        if os.path.exists(_path_to_lib):
            os.remove(_path_to_lib)
        raise NotImplementedError(f"The newest library code has not been compiled.")
    # Import the library.
    clib = ctypes.CDLL(_path_to_lib)
except:
    # Remove the shared object if it exists, because it is faulty.
    if os.path.exists(_shared_object_name):
        os.remove(_shared_object_name)
    # Compile a new shared object.
    _command = [_fort_compiler] + _ordered_dependencies + _compile_options + ["-o", _shared_object_name]
    if _verbose:
        print("Running system command with arguments")
        print("  ", " ".join(_command))
    # Run the compilation command.
    import subprocess
    subprocess.check_call(_command, cwd=_this_directory)
    # Import the shared object file as a C library with ctypes.
    clib = ctypes.CDLL(_path_to_lib)
# --------------------------------------------------------------------


# ----------------------------------------------
# Wrapper for the Fortran subroutine MATRIX_ROTATE

def matrix_rotate(thetax, thetay, thetaz):
    ''''''
    
    # Setting up "thetax"
    if (type(thetax) is not ctypes.c_float): thetax = ctypes.c_float(thetax)
    
    # Setting up "thetay"
    if (type(thetay) is not ctypes.c_float): thetay = ctypes.c_float(thetay)
    
    # Setting up "thetaz"
    if (type(thetaz) is not ctypes.c_float): thetaz = ctypes.c_float(thetaz)
    
    # Setting up "r"
    r = ctypes.c_void_p()
    r_dim_1 = ctypes.c_long()
    r_dim_2 = ctypes.c_long()

    # Call C-accessible Fortran wrapper.
    clib.c_matrix_rotate(ctypes.byref(thetax), ctypes.byref(thetay), ctypes.byref(thetaz), ctypes.byref(r_dim_1), ctypes.byref(r_dim_2), ctypes.byref(r))

    # Post-processing "r"
    r_size = (r_dim_1.value) * (r_dim_2.value)
    if (r_size > 0):
        r = numpy.array(ctypes.cast(r, ctypes.POINTER(ctypes.c_float*r_size)).contents, copy=False)
        r = r.reshape(r_dim_2.value,r_dim_1.value).T
    elif (r_size == 0):
        r = numpy.zeros(shape=(r_dim_2.value,r_dim_1.value), dtype=ctypes.c_float, order='F')
    else:
        r = None
    
    # Return final results, 'INTENT(OUT)' arguments only.
    return r


# ----------------------------------------------
# Wrapper for the Fortran subroutine MAIN_CALCULATION

def main_calculation(input_ints, a_list, theta_list):
    ''''''
    
    # Setting up "input_ints"
    if ((not issubclass(type(input_ints), numpy.ndarray)) or
        (not numpy.asarray(input_ints).flags.f_contiguous) or
        (not (input_ints.dtype == numpy.dtype(ctypes.c_int)))):
        import warnings
        warnings.warn("The provided argument 'input_ints' was not an f_contiguous NumPy array of type 'ctypes.c_int' (or equivalent). Automatically converting (probably creating a full copy).")
        input_ints = numpy.asarray(input_ints, dtype=ctypes.c_int, order='F')
    input_ints_dim_1 = ctypes.c_long(input_ints.shape[0])
    
    # Setting up "a_list"
    if ((not issubclass(type(a_list), numpy.ndarray)) or
        (not numpy.asarray(a_list).flags.f_contiguous) or
        (not (a_list.dtype == numpy.dtype(ctypes.c_float)))):
        import warnings
        warnings.warn("The provided argument 'a_list' was not an f_contiguous NumPy array of type 'ctypes.c_float' (or equivalent). Automatically converting (probably creating a full copy).")
        a_list = numpy.asarray(a_list, dtype=ctypes.c_float, order='F')
    a_list_dim_1 = ctypes.c_long(a_list.shape[0])
    a_list_dim_2 = ctypes.c_long(a_list.shape[1])
    
    # Setting up "theta_list"
    if ((not issubclass(type(theta_list), numpy.ndarray)) or
        (not numpy.asarray(theta_list).flags.f_contiguous) or
        (not (theta_list.dtype == numpy.dtype(ctypes.c_float)))):
        import warnings
        warnings.warn("The provided argument 'theta_list' was not an f_contiguous NumPy array of type 'ctypes.c_float' (or equivalent). Automatically converting (probably creating a full copy).")
        theta_list = numpy.asarray(theta_list, dtype=ctypes.c_float, order='F')
    theta_list_dim_1 = ctypes.c_long(theta_list.shape[0])
    theta_list_dim_2 = ctypes.c_long(theta_list.shape[1])

    # Call C-accessible Fortran wrapper.
    clib.c_main_calculation(ctypes.byref(input_ints_dim_1), ctypes.c_void_p(input_ints.ctypes.data), ctypes.byref(a_list_dim_1), ctypes.byref(a_list_dim_2), ctypes.c_void_p(a_list.ctypes.data), ctypes.byref(theta_list_dim_1), ctypes.byref(theta_list_dim_2), ctypes.c_void_p(theta_list.ctypes.data))

    # Return final results, 'INTENT(OUT)' arguments only.
    return 


# ----------------------------------------------
# Wrapper for the Fortran subroutine SENSOR_LINE_GENERATOR

def sensor_line_generator(x_boundary, nxf, pos_vector=None):
    ''''''
    
    # Setting up "x_boundary"
    if (type(x_boundary) is not ctypes.c_int): x_boundary = ctypes.c_int(x_boundary)
    
    # Setting up "nxf"
    if (type(nxf) is not ctypes.c_int): nxf = ctypes.c_int(nxf)
    
    # Setting up "pos_vector"
    if (pos_vector is None):
        pos_vector = numpy.zeros(shape=(3, nxf), dtype=ctypes.c_float, order='F')
    elif ((not issubclass(type(pos_vector), numpy.ndarray)) or
          (not numpy.asarray(pos_vector).flags.f_contiguous) or
          (not (pos_vector.dtype == numpy.dtype(ctypes.c_float)))):
        import warnings
        warnings.warn("The provided argument 'pos_vector' was not an f_contiguous NumPy array of type 'ctypes.c_float' (or equivalent). Automatically converting (probably creating a full copy).")
        pos_vector = numpy.asarray(pos_vector, dtype=ctypes.c_float, order='F')
    pos_vector_dim_1 = ctypes.c_long(pos_vector.shape[0])
    pos_vector_dim_2 = ctypes.c_long(pos_vector.shape[1])

    # Call C-accessible Fortran wrapper.
    clib.c_sensor_line_generator(ctypes.byref(x_boundary), ctypes.byref(nxf), ctypes.byref(pos_vector_dim_1), ctypes.byref(pos_vector_dim_2), ctypes.c_void_p(pos_vector.ctypes.data))

    # Return final results, 'INTENT(OUT)' arguments only.
    return pos_vector

