import numpy as np

def obj_fileTri2patch(file_in):
    """
    Reads a .obj triangular mesh file and outputs the vertex coordinates, faces, and related data.

    Parameters:
        file_in (str): Path to the .obj file.

    Returns:
        tuple: V (numpy.ndarray), F (numpy.ndarray), X (numpy.ndarray),
               Y (numpy.ndarray), Z (numpy.ndarray), M (numpy.ndarray)
    """
    vertices = []
    faces = []
    materials = []
    current_material = 0

    with open(file_in, 'r') as file:
        for line in file:
            if line.startswith('v '):
                # Vertex definition
                vertex = list(map(float, line.split()[1:]))
                vertices.append(vertex)
            elif line.startswith('usemtl '):
                # Material identifier
                current_material = int(line.split()[1].replace(';', ''))
            elif line.startswith('f '):
                # Face definition
                face_parts = line.split()[1:]
                face = []
                for part in face_parts:
                    # Handle different face definitions (v, v//vn, v/vt/vn)
                    indices = part.split('/')
                    face.append(int(indices[0]))
                faces.append(face)
                materials.append(current_material)

    # Convert to numpy arrays
    V = np.array(vertices)
    F = np.array(faces) - 1  # Convert to zero-based indexing
    M = np.array(materials)

    # Convert faces to X, Y, Z arrays
    X = V[F[:, 0], 0], V[F[:, 1], 0], V[F[:, 2], 0]
    Y = V[F[:, 0], 1], V[F[:, 1], 1], V[F[:, 2], 1]
    Z = V[F[:, 0], 2], V[F[:, 1], 2], V[F[:, 2], 2]

    return V, F, np.array(X), np.array(Y), np.array(Z), M
