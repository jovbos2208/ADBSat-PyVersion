import numpy as np
from .surfaceNormals import surface_normals
from .obj_fileTri2patch import obj_fileTri2patch
from scipy.io import savemat
import os

def ADBSatImport(file_in, path_out, struct_name, verbose=False):
    """
    Imports a triangular mesh from a .obj file, calculates relevant data,
    and saves the results into a .mat file.

    Parameters:
        file_in (str): Path to the input .obj file.
        path_out (str): Directory to save the output .mat file.
        struct_name (str): Name of the output structure.
        verbose (bool): If True, print progress information.

    Returns:
        str: Path to the saved .mat file.
    """
    if verbose:
        print(f"Processing file: {file_in}")

    # Read and process the .obj file
    vertices, faces, x_data, y_data, z_data, mat_id = obj_fileTri2patch(file_in)

    # Calculate surface normals, areas, and barycenters
    surface_normal, areas, barycenters = surface_normals(x_data, y_data, z_data)

    # Reference length
    len_ref = np.max(x_data) - np.min(x_data)

    # Create the mesh data structure
    meshdata = {
        'XData': x_data,
        'YData': y_data,
        'ZData': z_data,
        'MatID': mat_id,
        'Areas': areas,
        'SurfN': surface_normal,
        'BariC': barycenters,
        'Lref': len_ref
    }

    # Save to a .mat file
    file_out = os.path.join(path_out, f"{struct_name}.mat")
    savemat(file_out, {'meshdata': meshdata})

    if verbose:
        print(f"Mesh data saved to {file_out}")
        print("Summary:")
        print(f"  Number of elements: {len(areas)}")
        print(f"  Total area: {np.sum(areas):.6f}")
        print(f"  Maximum element area: {np.max(areas):.6f}")
        print(f"  Minimum element area: {np.min(areas):.6f}")
        print(f"  Reference length: {len_ref:.6f}")

    return file_out
