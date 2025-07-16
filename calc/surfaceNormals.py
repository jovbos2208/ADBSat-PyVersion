import numpy as np

def surface_normals(x, y, z):
    """
    Computes:
    - surface normals (unit vectors),
    - triangle areas,
    - barycenters (centroids) for a 3D triangle mesh.

    Parameters:
        x, y, z (3×N np.ndarray): Vertex coordinates per triangle

    Returns:
        tuple:
            - surfN: Normal vectors (3×N)
            - areas: Triangle areas (1×N)
            - bariC: Barycenters (3×N)
    """
    # Kantenvektoren
    V1 = np.array([x[1] - x[0], y[1] - y[0], z[1] - z[0]])
    V2 = np.array([x[2] - x[0], y[2] - y[0], z[2] - z[0]])

    # Normalen (ungewichtet)
    normals = np.cross(V1.T, V2.T).T
    magnitudes = np.linalg.norm(normals, axis=0)

    surfN = normals / magnitudes
    areas = 0.5 * magnitudes

    baryC = (np.array([x.sum(axis=0), y.sum(axis=0), z.sum(axis=0)])) / 3

    return surfN, areas, baryC

