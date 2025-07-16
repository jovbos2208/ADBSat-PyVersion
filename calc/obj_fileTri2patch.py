import numpy as np

def obj_fileTri2patch(file_in):
    """
    Parses a Wavefront .obj file into vertex and triangle arrays.

    Returns:
        tuple: vertices, faces, X, Y, Z coords (3xN), material IDs
    """
    vertices, faces, materials = [], [], []
    current_material = 0

    with open(file_in, 'r') as file:
        for line in file:
            if line.startswith('v '):
                vertices.append(list(map(float, line.split()[1:])))
            elif line.startswith('usemtl '):
                current_material = int(line.split()[1].replace(';', ''))
            elif line.startswith('f '):
                face = [int(part.split('/')[0]) for part in line.split()[1:]]
                faces.append(face)
                materials.append(current_material)

    V = np.array(vertices)
    F = np.array(faces) - 1  # OBJ uses 1-based indexing
    M = np.array(materials)

    X = V[F[:, 0], 0], V[F[:, 1], 0], V[F[:, 2], 0]
    Y = V[F[:, 0], 1], V[F[:, 1], 1], V[F[:, 2], 1]
    Z = V[F[:, 0], 2], V[F[:, 1], 2], V[F[:, 2], 2]

    return V, F, np.array(X), np.array(Y), np.array(Z), M

