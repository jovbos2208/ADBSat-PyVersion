import numpy as np

def insidetri(p1, p2, p3, points):
    """
    Prüft, ob Punkte innerhalb mehrerer Dreiecke liegen.

    Parameters:
        p1 (np.ndarray): Array mit den ersten Eckpunkten der Dreiecke, Form (2, N).
        p2 (np.ndarray): Array mit den zweiten Eckpunkten der Dreiecke, Form (2, N).
        p3 (np.ndarray): Array mit den dritten Eckpunkten der Dreiecke, Form (2, N).
        points (np.ndarray): Array der zu prüfenden Punkte, Form (2, M).

    Returns:
        np.ndarray: Boolean-Array mit True, wenn ein Punkt innerhalb eines Dreiecks liegt, Form (N, M).
    """
    def sign(p1, p2, p3):
        return (p1[0, :] - p3[0, :]) * (p2[1, :] - p3[1, :]) - (p2[0, :] - p3[0, :]) * (p1[1, :] - p3[1, :])

    # Anzahl der Dreiecke und Punkte
    num_triangles = p1.shape[1]
    num_points = points.shape[1]

    # Wiederhole Punkte, damit sie mit den Dreiecken verglichen werden können
    repeated_points = np.repeat(points[:, :, np.newaxis], num_triangles, axis=2)

    # Wiederhole Dreiecke, damit sie mit den Punkten verglichen werden können
    repeated_p1 = np.repeat(p1[:, :, np.newaxis], num_points, axis=2).transpose(0, 2, 1)
    repeated_p2 = np.repeat(p2[:, :, np.newaxis], num_points, axis=2).transpose(0, 2, 1)
    repeated_p3 = np.repeat(p3[:, :, np.newaxis], num_points, axis=2).transpose(0, 2, 1)

    # Berechne die Vorzeichen für die Punkte relativ zu den Dreieckskanten
    b1 = sign(repeated_points, repeated_p1, repeated_p2) < 0
    b2 = sign(repeated_points, repeated_p2, repeated_p3) < 0
    b3 = sign(repeated_points, repeated_p3, repeated_p1) < 0

    # Prüfe, ob alle Vorzeichen übereinstimmen (Punkt liegt im Dreieck)
    return np.logical_and.reduce([b1 == b2, b2 == b3], axis=0)




