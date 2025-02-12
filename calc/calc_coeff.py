import numpy as np
from .shadowAnaly import shadowAnaly
from .coeff_solar import coeff_solar
from .mainCoeff import mainCoeff
from scipy.io import loadmat, savemat
import os


def calc_coeff(fi_name, respath, aoaS, aosS, param_eq, flag_shad, flag_sol, delete_temp_files=False, verbose=False):
    """
    Calculates local and global coefficients for the triangular mesh geometry.

    Parameters:
        fi_name (str): Name of the .mat file containing the meshdata structure.
        respath (str): Path to save the results.
        aoaS (list): Angles of attack (radians).
        aosS (list): Angles of sideslip (radians).
        param_eq (dict): Dictionary of parameters for the selected GSI model.
        flag_shad (bool): Flag to perform shadow analysis.
        flag_sol (bool): Flag to calculate solar wind coefficients.
        delete_temp_files (bool): Flag to delete temporary files after merging.
        verbose (bool): If True, print progress.

    Returns:
        str: Path to the merged output file.
    """
    # Load mesh parameters
    mat_data = loadmat(fi_name)
    meshdata = mat_data['meshdata']
    x = meshdata['XData'][0, 0]
    y = meshdata['YData'][0, 0]
    z = meshdata['ZData'][0, 0]
    areas = meshdata['Areas'][0, 0]
    surfN = meshdata['SurfN'][0, 0]
    barC = meshdata['BariC'][0, 0]
    len_ref = meshdata['Lref'][0, 0]
    matID = meshdata['MatID'][0, 0]

    index_aoa = len(aoaS)
    index_aos = len(aosS)

    # Prepare output folder
    mat_name = os.path.splitext(os.path.basename(fi_name))[0]
    if index_aoa * index_aos > 1:
        folder_name = f"{mat_name}_{np.random.randint(1, 1000)}"
        os.makedirs(os.path.join(respath, folder_name), exist_ok=True)
        path_sav = os.path.join(respath, folder_name)
        is_batch = True
    else:
        path_sav = respath
        is_batch = False

    # Process each combination of aoa and aos
    for aoa in aoaS:
        for aos in aosS:

            # Calculate direction cosine matrices
            L_wb = np.array([
                [np.cos(aos) * np.cos(aoa), np.sin(aos), np.sin(aoa) * np.cos(aos)],
                [-np.sin(aos) * np.cos(aoa), np.cos(aos), -np.sin(aoa) * np.sin(aos)],
                [-np.sin(aoa), 0, np.cos(aoa)]
            ])

            L_gb = np.array([[1, 0, 0], [0, -1, 0], [0, 0, -1]])
            L_gw = np.dot(L_gb, L_wb.T)
            L_fb = np.array([[-1, 0, 0], [0, 1, 0], [0, 0, -1]])

            # Calculate angles between flow and normals
            vdir = np.dot(L_gw, np.array([-1, 0, 0]))
            vdir /= np.linalg.norm(vdir)
            v_matrix = np.tile(vdir[:, np.newaxis], surfN.shape[1])
            delta = np.arccos(np.einsum('ij,ij->j', -v_matrix, surfN))
            
            uD = v_matrix  # Einheitlicher Widerstandsvektor ist einfach die Strömungsrichtung
    
            # Berechnung des Lift-Vektors
            uL = -np.cross(np.cross(uD, surfN, axis=0), uD, axis=0)  # Kreuzprodukt berechnen
            uL_norm = np.linalg.norm(uL, axis=0)  # Norm berechnen
            
            # Vermeidung von Division durch Null (Panels normal zur Strömung)
            uL[:, uL_norm == 0] = -surfN[:, uL_norm == 0]  # Setze Lift-Vektor gleich Oberflächennormalen, wenn undefiniert
            valid_idx = uL_norm > 0  # Indizes der gültigen Vektoren
            
            # Normiere gültige Lift-Vektoren
            uL[:, valid_idx] /= uL_norm[valid_idx]

            # Berechnung von `gamma` und `ell`
            param_eq["gamma"] = np.einsum('ij,ij->j', -uD, surfN)  # Elementweises Skalarprodukt über Spalten
            param_eq["ell"] = np.einsum('ij,ij->j', -uL, surfN)  # Dasselbe für uL


            # Calculate aerodynamic coefficients
            cp, ctau, cd, cl = mainCoeff(param_eq, delta, matID)

            # Optional: Solar coefficients
            if flag_sol:
                cn, cs = coeff_solar(delta, param_eq)

            # Shadow analysis
            if flag_shad:
                shad_pan = shadowAnaly(x, y, z, barC, delta, L_gw)
                cp[shad_pan] = ctau[shad_pan] = cd[shad_pan] = cl[shad_pan] = 0

            # Save results
            method_name = param_eq['gsi_model']
            file_name = f"{mat_name}_{method_name}_aoa{int(np.degrees(aoa))}_aos{int(np.degrees(aos))}.mat"
            savemat(os.path.join(path_sav, file_name), {
                'cp': cp, 'ctau': ctau, 'cd': cd, 'cl': cl,
                'aoa': aoa, 'aos': aos, 'delta': delta
            })

            path_sav = os.path.join(path_sav, file_name)

            if verbose:
                print(f"Saved results for AoA: {np.degrees(aoa)} deg, AoS: {np.degrees(aos)} deg")

    return path_sav
