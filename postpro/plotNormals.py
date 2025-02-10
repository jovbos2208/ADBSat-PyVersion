import numpy as np
import scipy.io
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def plot_normals(fi_name):
    """
    Plots the mesh and normal vectors.
    
    Parameters:
        fi_name (str): Path to the .mat file containing the meshdata structure.
    """
    # Load .mat file
    data = scipy.io.loadmat(fi_name)
    meshdata = data['meshdata']
    
    # Extract mesh data
    x = meshdata['XData'][0, 0]
    y = meshdata['YData'][0, 0]
    z = meshdata['ZData'][0, 0]
    barC = meshdata['BariC'][0, 0]
    surfN = meshdata['SurfN'][0, 0]
    matID = meshdata['MatID'][0, 0].flatten()
    
    # Generate colormap
    unique_mats = np.unique(matID)
    mats = len(unique_mats)
    cmap = plt.cm.get_cmap('tab10', mats)
    
    # First figure: Normals plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.quiver(barC[0, :], barC[1, :], barC[2, :], 
              surfN[0, :], surfN[1, :], surfN[2, :], length=0.1, normalize=True)
    
    # Mesh plot
    for i in range(len(x)):
        verts = [list(zip(x[i], y[i], z[i]))]
        poly = Poly3DCollection(verts, alpha=0.5, edgecolor='k')
        ax.add_collection3d(poly)
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('ADBSat Mesh')
    ax.set_box_aspect([1,1,1])
    plt.show()
    
    # Second figure: Material ID plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    for i in range(len(x)):
        verts = [list(zip(x[i], y[i], z[i]))]
        poly = Poly3DCollection(verts, alpha=0.7, edgecolor='none')
        poly.set_facecolor(cmap(matID[i] / mats))
        ax.add_collection3d(poly)
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('ADBSat Material ID')
    
    # Colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0, vmax=mats))
    cb = fig.colorbar(sm, ax=ax, ticks=np.arange(mats + 1))
    cb.set_label("Material ID")
    
    ax.set_box_aspect([1,1,1])
    plt.show()


