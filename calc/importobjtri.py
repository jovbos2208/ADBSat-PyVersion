import os
import numpy as np
import scipy.io
from .surfaceNormals import surface_normals
from .obj_fileTri2patch import obj_fileTri2patch


def importobjtri(file_in, path_out, stru_name, verbose=False):
    """
    Imports a triangular surface mesh from a Wavefront OBJ file, computes geometric properties,
    and saves the result in a MATLAB .mat structure for further processing.

    Parameters:
        file_in (str): Path to the input .obj file (triangular mesh format).
        path_out (str): Directory where the output .mat file will be saved.
        stru_name (str): Name of the output structure (used in the file name).
        verbose (bool): If True, print import statistics and file info.

    Returns:
        str: Full path to the saved .mat file.
    """
    # Step 1: Load vertices and triangles from OBJ file
    vertices, faces, X, Y, Z, materials = obj_fileTri2patch(file_in)

    # Step 2: Compute surface normals, panel areas, and barycenters
    surface_normal, areas, barycenters = surface_normals(X, Y, Z)

    # Step 3: Calculate reference length (bounding box in X-direction)
    Lref = np.max(X) - np.min(X)

    # Step 4: Create mesh structure to save
    meshdata = {
        'XData': X,
        'YData': Y,
        'ZData': Z,
        'MatID': materials,
        'Areas': areas,
        'SurfN': surface_normal,
        'BariC': barycenters,
        'Lref': Lref
    }

    # Step 5: Save mesh structure to MATLAB file
    output_path = os.path.join(path_out, f"{stru_name}.mat")
    scipy.io.savemat(output_path, {'meshdata': meshdata})

    # Optional: print detailed statistics
    if verbose:
        num_faces = X.shape[1]
        total_area = np.sum(areas)
        max_area = np.max(areas)
        min_area = np.min(areas)
        num_materials = np.max(materials)

        print("✅ Import completed!")
        print("======================================")
        print(f"Triangles:              {num_faces}")
        print(f"Total surface area:     {total_area:.6f}")
        print(f"Maximum element area:   {max_area:.6f}")
        print(f"Minimum element area:   {min_area:.6f}")
        print(f"Reference length (X):   {Lref:.6f}")
        print(f"Material ID range:      {num_materials}")
        print(f"Saved as:               {output_path}")

    return output_path

