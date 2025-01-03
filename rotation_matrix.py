import numpy as np

def rotation_total(thetax, thetay, thetaz):
    # Implement inline matrix multiplication for performance.
    return np.array([
        [np.cos(thetay) * np.cos(thetaz), -np.cos(thetay) * np.sin(thetaz), np.sin(thetay)],
        [np.cos(thetax) * np.sin(thetaz) + np.sin(thetax) * np.sin(thetay) * np.cos(thetaz),
         np.cos(thetax) * np.cos(thetaz) - np.sin(thetax) * np.sin(thetay) * np.sin(thetaz),
         -np.sin(thetax) * np.cos(thetay)],
        [np.sin(thetax) * np.sin(thetaz) - np.cos(thetax) * np.sin(thetay) * np.cos(thetaz),
         np.sin(thetax) * np.cos(thetaz) + np.cos(thetax) * np.sin(thetay) * np.sin(thetaz),
         np.cos(thetax) * np.cos(thetay)]
    ])
    