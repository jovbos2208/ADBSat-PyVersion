import numpy as np
from .surfaceNormals import surface_normals
from .obj_fileTri2patch import obj_fileTri2patch
from scipy.io import savemat
import os


def ADBSatImport(file_in, path_out, struct_name, verbose=False):
    """
    Imports a triangular surface mesh from a Wavefront OBJ file,
    calculates normals, areas, barycenters, and writes the data to a MATLAB .mat file.

    Parameters:
        file_in (str): Path to the input .obj file.
        path_out (str): Directory for saving the .mat output.
        struct_name (str): Variable name for the output structure (in the .mat file).
        verbose (bool): If True, print detailed processing info.

    Returns:
        str: Path to the saved .mat file.
    """
    if verbose:
        print(f"Processing mesh: {file_in}")

    # 1. Read OBJ file into vertex & triangle patch representation
    vertices, faces, x_data, y_data, z_data, mat_id = obj_fileTri2patch(file_in)

    # 2. Compute geometric properties
    surface_normal, areas, barycenters = surface_normals(x_data, y_data, z_data)

    # 3. Reference length (used for normalization and scaling)
    len_ref = np.max(x_data) - np.min(x_data)

    # 4. Create output dictionary
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

    # 5. Save to .mat file
    file_out = os.path.join(path_out, f"{struct_name}.mat")
    savemat(file_out, {'meshdata': meshdata})

    if verbose:
        print(f"Mesh saved to: {file_out}")
        print("Mesh summary:")
        print(f"  Triangles:           {len(areas)}")
        print(f"  Total area:          {np.sum(areas):.6f}")
        print(f"  Largest element:     {np.max(areas):.6f}")
        print(f"  Smallest element:    {np.min(areas):.6f}")
        print(f"  Reference length:    {len_ref:.6f}")

    return file_out

