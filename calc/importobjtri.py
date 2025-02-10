import os
import numpy as np
import scipy.io
from .surfaceNormals import surface_normals
from .obj_fileTri2patch import obj_fileTri2patch

def importobjtri(file_in, path_out, stru_name, verbose=False):
    """
    Imports a triangular mesh from a .obj file and saves it in .mat format.
    
    Parameters:
        file_in (str): Input filepath of the .obj file.
        path_out (str): Directory to save the output .mat file.
        stru_name (str): Output name for the model file (without extension).
        verbose (bool): Flag for printing detailed output to the console.
    
    Returns:
        str: Full path to the saved .mat file.
    """
    # Load the triangular mesh
    vertices, faces, X, Y, Z, materials = obj_fileTri2patch(file_in)
    
    # Calculate surface normals, areas, and barycenters
    surface_normals, areas, barycenters = surface_normals(X, Y, Z)
    
    # Compute reference length
    Lref = np.max(X) - np.min(X)
    
    # Create the mesh data dictionary
    meshdata = {
        'XData': X,
        'YData': Y,
        'ZData': Z,
        'MatID': materials,
        'Areas': areas,
        'SurfN': surface_normals,
        'BariC': barycenters,
        'Lref': Lref
    }
    
    # Save the data in .mat format
    output_path = os.path.join(path_out, f"{stru_name}.mat")
    scipy.io.savemat(output_path, {'meshdata': meshdata})
    
    # Extract and print mesh statistics if verbose mode is on
    num_faces = X.shape[1]
    total_area = np.sum(areas)
    max_area = np.max(areas)
    min_area = np.min(areas)
    num_materials = np.max(materials)
    
    if verbose:
        print("Import finished!")
        print("******************************************")
        print(f"Number of elements: {num_faces}")
        print(f"Total area: {total_area}")
        print(f"Maximum element area: {max_area}")
        print(f"Minimum element area: {min_area}")
        print(f"Reference length (maxX-minX): {Lref}")
        print(f"Number of material references: {num_materials}")
        print(f"Created file '{output_path}'")
    
    return output_path
