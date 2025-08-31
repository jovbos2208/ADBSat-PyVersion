import numpy as np
from .shadowAnaly import shadowAnaly
from .coeff_solar import coeff_solar
from .mainCoeff import mainCoeff
from scipy.io import loadmat, savemat
import os


def calc_coeff(fi_name, respath, aoaS, aosS, param_eq,
               flag_shad, flag_sol, dp,
               delete_temp_files=False, verbose=False):
    """
    Calculates local and global aerodynamic and (optional) solar radiation pressure coefficients
    for a triangular surface mesh using the specified GSI model.

    Parameters:
        fi_name (str): Path to the input .mat file with meshdata structure.
        respath (str): Output folder for result files.
        aoaS (list[float]): List of angles of attack [rad].
        aosS (list[float]): List of angles of sideslip [rad].
        param_eq (dict): GSI model parameters (Tw, Tinf, vinf, etc.).
        flag_shad (bool): Enable shadowing analysis (zero contribution from occluded panels).
        flag_sol (bool): Enable solar radiation pressure coefficient computation.
        dp: Placeholder for compatibility.
        delete_temp_files (bool): Currently unused; reserved for cleanup logic.
        verbose (bool): Print progress logs.

    Returns:
        str: Path to the last generated .mat file.
    """

    # --------------------------
    # Load mesh data
    # --------------------------
    mat_data = loadmat(fi_name)
    meshdata = mat_data['meshdata']

    x = meshdata['XData'][0, 0]
    y = meshdata['YData'][0, 0]
    z = meshdata['ZData'][0, 0]
    areas = meshdata['Areas'][0, 0].ravel()
    surfN = meshdata['SurfN'][0, 0]
    barC = meshdata['BariC'][0, 0]
    len_ref = float(meshdata['Lref'][0, 0])
    matID = meshdata['MatID'][0, 0].ravel()

    def _repeat(v, n):
        return np.tile(v.reshape(3, 1), n)

    # Output folder setup
    idx_aoa = len(aoaS)
    idx_aos = len(aosS)
    mat_name = os.path.splitext(os.path.basename(fi_name))[0]

    if idx_aoa * idx_aos > 1:
        folder_name = f"{mat_name}_{np.random.randint(1, 1000):03d}"
        os.makedirs(os.path.join(respath, folder_name), exist_ok=True)
        path_res = os.path.join(respath, folder_name)
    else:
        path_res = respath

    # --------------------------
    # Loop over angle combinations
    # --------------------------
    for aoa in aoaS:
        for aos in aosS:
            # --- Direction cosine matrices ---
            L_wb = np.array([
                [np.cos(aos) * np.cos(aoa),  np.sin(aos),              np.sin(aoa) * np.cos(aos)],
                [-np.sin(aos) * np.cos(aoa), np.cos(aos),             -np.sin(aoa) * np.sin(aos)],
                [-np.sin(aoa),               0.0,                      np.cos(aoa)]
            ])
            L_gb = np.diag([1.0, -1.0, -1.0])  # ground → body
            L_gw = L_gb @ L_wb.T              # ground → wind
            L_fb = np.diag([-1.0, 1.0, -1.0])  # flow → body

            # --- Flow direction in wind frame ---
            v_dir = L_gw @ np.array([-1.0, 0.0, 0.0])
            v_dir /= np.linalg.norm(v_dir)
            v_matrix = _repeat(v_dir, surfN.shape[1])

            # --- Angle between surface normal and freestream direction ---
            delta = np.arccos(np.einsum('ij,ij->j', -v_matrix, surfN))

            # --- Update model parameters with local angles ---
            param_eq = param_eq.copy()
            param_eq["gamma"] = np.cos(delta)
            param_eq["ell"] = np.sin(delta)

            # --- Local aerodynamic coefficients (per panel) ---
            cp, ctau, cd, cl = mainCoeff(param_eq, delta, matID)

            # --- Optional solar pressure model ---
            if flag_sol:
                cn, cs = coeff_solar(delta, param_eq)

            # --- Shadowing analysis ---
            if flag_shad:
                shad_pan = shadowAnaly(x, y, z, barC, delta, L_gw)
                cp[shad_pan] = ctau[shad_pan] = cd[shad_pan] = cl[shad_pan] = 0.0
                if flag_sol:
                    cn[shad_pan] = cs[shad_pan] = 0.0

            # ============================
            # GLOBAL FORCES & MOMENTS
            # ============================

            # Reference areas
            area_proj = areas * np.cos(delta)
            area_total = np.sum(areas)
            area_ref = area_total / 2.0

            # Shear direction unit vectors
            tau_dir = np.cross(surfN.T, np.cross(v_matrix.T, surfN.T)).T
            tau_norm = np.linalg.norm(tau_dir, axis=0)
            good = tau_norm > 0
            tau_dir[:, good] /= tau_norm[good]
            tau_dir[:, ~good] = 0.0

            # --- Global aerodynamic force ---
            ctau_area = ctau * areas
            cp_area = cp * areas
            Cf_g = (tau_dir @ ctau_area - surfN @ cp_area) / area_ref
            Cf_w = L_gw.T @ Cf_g
            Cf_f = L_fb @ (L_gb.T @ Cf_g)

            # --- Global aerodynamic moment ---
            cross_tau = np.cross(barC.T, tau_dir.T).T
            cross_nrm = np.cross(barC.T, -surfN.T).T
            Cm_g = (cross_tau @ ctau_area + cross_nrm @ cp_area) / (area_ref * len_ref)
            Cm_b = L_gb.T @ Cm_g

            print(f'AoS: {np.degrees(aos)-90:5.1f}, Cd-Value: {np.mean(cd):5.5f}')

            # --- Global solar coefficients (if enabled) ---
            if flag_sol:
                cstau = cs * np.sin(delta)
                csn = cs * np.cos(delta)
                cstau_area = cstau * areas
                csn_area = (cn + csn) * areas

                Cs_g = (tau_dir @ cstau_area - surfN @ csn_area) / area_ref
                Cf_s = L_gw.T @ Cs_g

                Cm_s_g = (cross_tau @ cstau_area + cross_nrm @ csn_area) / (area_ref * len_ref)
                Cm_s_b = L_gb.T @ Cm_s_g
            else:
                Cs_g = Cf_s = Cm_s_g = Cm_s_b = np.zeros(3)

            # --------------------------
            # Save to .mat file
            # --------------------------
            method = param_eq['gsi_model']
            fn_out = f"{mat_name}_{method}_aoa{int(np.degrees(aoa))}_aos{int(np.degrees(aos))}.mat"
            path_out = os.path.join(path_res, fn_out)

            savemat(path_out, {
                'cp': cp, 'ctau': ctau, 'cd': cd, 'cl': cl,
                'delta': delta,
                'Cf_g': Cf_g, 'Cf_w': Cf_w, 'Cf_f': Cf_f,
                'Cm_G': Cm_g, 'Cm_B': Cm_b,
                'Cs_G': Cs_g, 'Cf_s': Cf_s,
                'Cm_s_G': Cm_s_g, 'Cm_s_B': Cm_s_b,
                'aoa': aoa, 'aos': aos,
                'AreaRef': area_ref,
                'AreaTotal': area_total,
                'AreaProj': area_proj,
                'LenRef': len_ref
            })

            if verbose:
                print(f"[calc_coeff] AoA = {np.degrees(aoa):5.1f}°, AoS = {np.degrees(aos):5.1f}° → {fn_out}")

            last_path = path_out

    # Optional: delete temp files
    if delete_temp_files:
        pass  # placeholder

    return last_path

