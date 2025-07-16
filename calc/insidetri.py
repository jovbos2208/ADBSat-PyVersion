import numpy as np

def insidetri(p1, p2, p3, points):
    """
    Tests whether multiple 2D points lie inside multiple triangles.

    Parameters:
        p1, p2, p3 (np.ndarray): Triangle corners, shape (2, N)
        points (np.ndarray): Points to test, shape (2, M)

    Returns:
        np.ndarray: Boolean array (N, M), True where point lies inside triangle
    """
    def sign(p1, p2, p3):
        return (p1[0, :] - p3[0, :]) * (p2[1, :] - p3[1, :]) - (p2[0, :] - p3[0, :]) * (p1[1, :] - p3[1, :])

    num_triangles = p1.shape[1]
    num_points = points.shape[1]

    # Broadcast triangles and points
    points_rep = np.repeat(points[:, :, np.newaxis], num_triangles, axis=2)
    p1_rep = np.repeat(p1[:, :, np.newaxis], num_points, axis=2).transpose(0, 2, 1)
    p2_rep = np.repeat(p2[:, :, np.newaxis], num_points, axis=2).transpose(0, 2, 1)
    p3_rep = np.repeat(p3[:, :, np.newaxis], num_points, axis=2).transpose(0, 2, 1)

    # Sign tests
    b1 = sign(points_rep, p1_rep, p2_rep) < 0
    b2 = sign(points_rep, p2_rep, p3_rep) < 0
    b3 = sign(points_rep, p3_rep, p1_rep) < 0

    return np.logical_and.reduce([b1 == b2, b2 == b3], axis=0)

