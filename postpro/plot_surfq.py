import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import mplcursors

def plot_surfq(file_in, mod_in, aoa_deg, aos_deg, param, save_path=None):
    """
    Plots the surface mesh with color proportional to the chosen parameter.

    Parameters:
        file_in (str): Path to the file containing the results.
        mod_in (str): Path to the file containing the mesh data.
        aoa_deg (float): Angle of attack in degrees.
        aos_deg (float): Angle of sideslip in degrees.
        param (str): Surface parameter to plot (e.g., 'cp', 'ctau', 'cd', 'cl').
        save_path (str, optional): Path to save the plotted figure. If None, the figure is not saved.

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
    for i in range(num_faces):
        verts.append([(x[j, i], y[j, i], z[j, i]) for j in range(3)])

    # Create the plot
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Add mesh with color-coded parameter values
    collection = Poly3DCollection(verts, cmap='viridis', edgecolor='k', linewidth=0.5)
    collection.set_array(param_values)
    ax.add_collection3d(collection)

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
