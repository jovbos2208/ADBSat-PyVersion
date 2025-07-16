import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import mplcursors

def plot_surfq(file_in, mod_in, aoa_deg, aos_deg, param, save_path=None, show_normals=True, normal_scale=0.1):
    """
    Color-coded surface visualization for a chosen aerodynamic parameter.

    Parameters:
        file_in (str): Path to result .mat file.
        mod_in (str): Path to mesh .mat file.
        aoa_deg (float): Angle of attack [deg]
        aos_deg (float): Angle of sideslip [deg]
        param (str): Parameter to plot ('cp', 'ctau', 'cd', etc.)
        save_path (str): Optional path to save the plot.
        show_normals (bool): Whether to show normals.
        normal_scale (float): Scaling for normal vectors.
    """
    mesh_data = loadmat(mod_in)['meshdata']
    x, y, z = mesh_data['XData'][0, 0], mesh_data['YData'][0, 0], mesh_data['ZData'][0, 0]
    results = loadmat(file_in)

    if param not in results:
        raise KeyError(f"'{param}' not found in result file.")

    param_values = results[param].flatten()

    verts, centers, normals = [], [], []
    for i in range(x.shape[1]):
        v0 = np.array([x[0, i], y[0, i], z[0, i]])
        v1 = np.array([x[1, i], y[1, i], z[1, i]])
        v2 = np.array([x[2, i], y[2, i], z[2, i]])
        verts.append([v0, v1, v2])
        center = (v0 + v1 + v2) / 3.0
        centers.append(center)
        n = np.cross(v1 - v0, v2 - v0)
        normals.append(n / np.linalg.norm(n) if np.linalg.norm(n) else n)

    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection='3d')
    collection = Poly3DCollection(verts, cmap='viridis', edgecolor='k', linewidth=0.5)
    collection.set_array(param_values)
    ax.add_collection3d(collection)

    if show_normals:
        c = np.array(centers)
        n = np.array(normals)
        ax.quiver(c[:, 0], c[:, 1], c[:, 2], n[:, 0], n[:, 1], n[:, 2],
                  length=normal_scale, color='red', linewidth=0.8)

    ax.set_title(f"{param.upper()} Distribution (AoA={aoa_deg}°, AoS={aos_deg}°)")
    plt.colorbar(collection, ax=ax, pad=0.1, label=f"{param.upper()} Coefficient")

    mplcursors.cursor(collection, hover=True).connect(
        "add", lambda sel: sel.annotation.set_text(f"{param}: {param_values[sel.index]:.3f}")
    )

    if save_path:
        plt.savefig(save_path, dpi=300)
        print(f"Saved to {save_path}")
    plt.show()

