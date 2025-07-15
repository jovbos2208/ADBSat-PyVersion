import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import mplcursors

def plot_surfq(file_in, mod_in, aoa_deg, aos_deg, param, save_path=None, show_normals=True, normal_scale=0.1):
    """
    Plots the surface mesh with color proportional to the chosen parameter and optional surface normals.

    Parameters:
        file_in (str): Path to the file containing the results.
        mod_in (str): Path to the file containing the mesh data.
        aoa_deg (float): Angle of attack in degrees.
        aos_deg (float): Angle of sideslip in degrees.
        param (str): Surface parameter to plot (e.g., 'cp', 'ctau', 'cd', 'cl').
        save_path (str, optional): Path to save the plotted figure. If None, the figure is not saved.
        show_normals (bool): Whether to show surface normals.
        normal_scale (float): Scale factor for normal vector length.

    Returns:
        None
    """
    # Load model mesh data
    mesh_data = loadmat(mod_in)['meshdata']
    x = mesh_data['XData'][0, 0]
    y = mesh_data['YData'][0, 0]
    z = mesh_data['ZData'][0, 0]

    # Load aerodynamic results
    results = loadmat(file_in)
    if 'aedb' in results:
        raise ValueError("Please select a single ADBSat output .mat file.")

    # Extract parameter values
    if param not in results:
        raise KeyError(f"Parameter '{param}' not found in the results file.")
    param_values = results[param].flatten()

    # Prepare vertex coordinates for the mesh
    num_faces = x.shape[1]
    verts = []
    face_centers = []
    normals = []

    for i in range(num_faces):
        v0 = np.array([x[0, i], y[0, i], z[0, i]])
        v1 = np.array([x[1, i], y[1, i], z[1, i]])
        v2 = np.array([x[2, i], y[2, i], z[2, i]])
        verts.append([v0, v1, v2])

        # Compute normal and center for each face
        edge1 = v1 - v0
        edge2 = v2 - v0
        normal = np.cross(edge1, edge2)
        norm_len = np.linalg.norm(normal)
        if norm_len > 0:
            normal /= norm_len
        center = (v0 + v1 + v2) / 3.0
        face_centers.append(center)
        normals.append(normal)

    # Convert to numpy arrays
    face_centers = np.array(face_centers)
    normals = np.array(normals)

    # Create the plot
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Add mesh with color-coded parameter values
    collection = Poly3DCollection(verts, cmap='viridis', edgecolor='k', linewidth=0.5)
    collection.set_array(param_values)
    ax.add_collection3d(collection)

    # Show normals as arrows
    if show_normals:
        ax.quiver(
            face_centers[:, 0], face_centers[:, 1], face_centers[:, 2],
            normals[:, 0], normals[:, 1], normals[:, 2],
            length=normal_scale, color='red', linewidth=0.8, normalize=True
        )

    # Set axis limits
    ax.set_xlim([x.min(), x.max()])
    ax.set_ylim([y.min(), y.max()])
    ax.set_zlim([z.min(), z.max()])

    # Add color bar
    cbar = plt.colorbar(collection, ax=ax, pad=0.1)
    cbar.set_label(f"{param} Coefficients")

    # Set labels and title
    ax.set_title(f"{param} Surface Distribution\nAoA: {aoa_deg:.2f} deg, AoS: {aos_deg:.2f} deg")
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    # Add interactivity with mplcursors
    mplcursors.cursor(collection, hover=True).connect(
        "add", lambda sel: sel.annotation.set_text(f"{param}: {param_values[sel.index]:.3f}")
    )

    # Save the figure if save_path is provided
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Figure saved to {save_path}")

    # Show the plot
    plt.show()

