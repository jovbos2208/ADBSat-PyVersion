import numpy as np

def surface_normals(x, y, z):
    """
    Calculates the surface normals, areas, and barycenters of a triangular mesh.

    Parameters:
        x (ndarray): X coordinates of the vertices of the mesh elements (3xN).
        y (ndarray): Y coordinates of the vertices of the mesh elements (3xN).
        z (ndarray): Z coordinates of the vertices of the mesh elements (3xN).

    Returns:
        tuple: surfN (3xN normalized normal vectors),
               areas (1xN areas of each triangle),
               bariC (3xN barycenters of each triangle).
    """
    # Calculate vectors for each triangle
    xV1 = x[1, :] - x[0, :]
    xV2 = x[2, :] - x[0, :]

    yV1 = y[1, :] - y[0, :]
    yV2 = y[2, :] - y[0, :]

    zV1 = z[1, :] - z[0, :]
    zV2 = z[2, :] - z[0, :]

    V1 = np.array([xV1, yV1, zV1])
    V2 = np.array([xV2, yV2, zV2])

    # Calculate barycenters
    bariC = (np.array([x[0, :] + x[1, :] + x[2, :], 
                       y[0, :] + y[1, :] + y[2, :], 
                       z[0, :] + z[1, :] + z[2, :]]) / 3)

    # Calculate surface normals
    surfN = np.cross(V1.T, V2.T).T
    ModN = np.linalg.norm(surfN, axis=0)
    surfN = surfN / ModN

    # Calculate areas
    areas = 0.5 * ModN

    return surfN, areas, bariC
