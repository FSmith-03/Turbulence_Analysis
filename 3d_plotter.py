import numpy as np
import plotly.graph_objects as go
import rotation_matrix as k

def plot3d(N_E):
    L = 1.0
    x_boundary = 3
    y_boundary = 3
    z_boundary = 10
    Nx = 50
    Nyz = 50

    x, y, z = np.meshgrid(np.linspace(-x_boundary*L, x_boundary*L, Nx), 
                          np.linspace(-y_boundary*L, y_boundary*L, Nyz), 
                          np.linspace(-z_boundary*L, z_boundary*L, Nyz))
    u = np.zeros_like(x)
    v = np.zeros_like(y)
    w = np.zeros_like(z)
    # Flatten the meshgrid arrays
    points = np.vstack([x.ravel(), y.ravel(), z.ravel()])
    for i in range(N_E):
        a = np.array([0, 0, 0])
        theta_x, theta_y, theta_z = np.array([0, 0, 0])
        # Combine rotation matrices
        R = k.rotation_total(theta_x, theta_y, theta_z)
        
        # Apply the combined rotation matrix
        rotated_points = R @ points
        
        # Translate the rotated points by vector a
        translated_points = rotated_points + np.array(a).reshape(3, 1)
        
        # Reshape back to the original shape
        x_r = translated_points[0].reshape(x.shape)
        y_r = translated_points[1].reshape(y.shape)
        z_r = translated_points[2].reshape(z.shape)
        
        #u_0 = -y_r * np.exp(-(x_r**2 + y_r**2 + z_r**2) / (2*L**2))
        #v_0 = x_r * np.exp(-(x_r**2 + y_r**2 + z_r**2) / (2*L**2))
        #w_0 = np.zeros_like(u)
        ARN = 1
        ARP = 1
        u_0 = -y_r * np.exp(-2 * (((x_r**2 + y_r**2) / (ARP * L)**2) + (z_r**2 / (ARN * L)**2)))
        v_0 = x_r * np.exp(-2 * (((x_r**2 + y_r**2) / (ARP * L)**2) + (z_r**2 / (ARN * L)**2)))
        w_0 = np.zeros_like(u_0)
        u = np.add(u, u_0)
        v = np.add(v, v_0)
        w = np.add(w, w_0)
    mag_v = np.sqrt(u**2 + v**2 + w**2)
    # Create isosurface plot
    fig = go.Figure(data=go.Isosurface(
        x=x.flatten(),
        y=y.flatten(),
        z=z.flatten(),
        value=mag_v.flatten(),
        isomin=mag_v.min(),
        isomax=mag_v.max(),
        surface_count=10,  # Number of isosurfaces
        caps=dict(x_show=False, y_show=False, z_show=False)
    ))
    fig.update_layout(scene=dict(
        xaxis_title='X',
        yaxis_title='Y',
        zaxis_title='Z'
    ))
    
    fig.show()

plot3d(1)