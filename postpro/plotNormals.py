import numpy as np
import scipy.io
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def plot_normals(fi_name):
    """
    Visualizes the mesh and surface normals from a .mat file.

    Parameters:
        fi_name (str): Path to the .mat file containing the 'meshdata' structure.
    """
    data = scipy.io.loadmat(fi_name)
    meshdata = data['meshdata']

    x, y, z = meshdata['XData'][0, 0], meshdata['YData'][0, 0], meshdata['ZData'][0, 0]
    barC = meshdata['BariC'][0, 0]
    surfN = meshdata['SurfN'][0, 0]
    matID = meshdata['MatID'][0, 0].flatten()

    unique_mats = np.unique(matID)
    cmap = plt.cm.get_cmap('tab10', len(unique_mats))

    # Normals plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.quiver(barC[0], barC[1], barC[2], surfN[0], surfN[1], surfN[2], length=0.1)

    for i in range(len(x)):
        verts = [list(zip(x[i], y[i], z[i]))]
        poly = Poly3DCollection(verts, alpha=0.5, edgecolor='k')
        ax.add_collection3d(poly)

    ax.set_title('ADBSat Mesh with Normals')
    plt.show()

    # Material ID plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for i in range(len(x)):
        verts = [list(zip(x[i], y[i], z[i]))]
        poly = Poly3DCollection(verts, alpha=0.7)
        poly.set_facecolor(cmap(matID[i] / len(unique_mats)))
        ax.add_collection3d(poly)

    ax.set_title('Mesh Colored by Material ID')
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0, vmax=len(unique_mats)))
    fig.colorbar(sm, ax=ax, ticks=np.arange(len(unique_mats) + 1), label='Material ID')
    plt.show()

