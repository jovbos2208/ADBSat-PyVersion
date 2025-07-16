import os
import numpy as np
from math import radians
from .ADBSatConstants import EnvironmentData, ConstantsData
from .environment import environment
from .calc_coeff import calc_coeff


def adbsat_fcn(mod_path, res_path, param_eq, aoa_deg, aos_deg,
               flag_shadow, flag_solar, env, del_files, verbose):
    """
    Main interface function for aerodynamic & solar coefficient analysis.

    This function prepares inputs, validates parameters, and orchestrates the computation
    of local/global aerodynamic and solar radiation pressure (SRP) coefficients based on
    the selected GSI model.

    Parameters:
        mod_path (str): Path to the input mesh model (*.mat file).
        res_path (str): Output directory for result files.
        param_eq (dict): Parameter dictionary for GSI and SRP models.
        aoa_deg (list[float]): List of angles of attack (in degrees).
        aos_deg (list[float]): List of angles of sideslip (in degrees).
        flag_shadow (bool): If True, apply shadow analysis.
        flag_solar (bool): If True, compute solar pressure coefficients.
        env (list): Environmental context (e.g., MSISE values).
        del_files (bool): If True, delete temporary files (currently not used).
        verbose (bool): If True, print additional information.

    Returns:
        str: Path to the last generated result file (*.mat).
    """
    # ----------------------------
    # Validate GSI model inputs
    # ----------------------------
    gsi_model = param_eq.get('gsi_model', '').lower()

    if gsi_model == 'cll':
        if not all(key in param_eq for key in ['alphaN', 'sigmaT']):
            raise ValueError("CLL model requires 'alphaN' and 'sigmaT'.")

    elif gsi_model in ['schaaf', 'storchhyp']:
        if not all(key in param_eq for key in ['sigmaN', 'sigmaT']):
            raise ValueError(f"{gsi_model} model requires 'sigmaN' and 'sigmaT'.")

    elif gsi_model in ['cook', 'sentman', 'maxwell', 'dria']:
        if 'alpha' not in param_eq:
            raise ValueError(f"{gsi_model} model requires 'alpha'.")

    elif gsi_model == 'newton':
        pass  # No parameters needed

    else:
        raise ValueError("Unrecognized GSI model name in 'param_eq'.")

    # ----------------------------
    # Validate solar inputs
    # ----------------------------
    if flag_solar and not all(key in param_eq for key in ['sol_cR', 'sol_cD']):
        raise ValueError("Solar model requires 'sol_cR' and 'sol_cD' in 'param_eq'.")

    # ----------------------------
    # Default values and pre-processing
    # ----------------------------
    if 'Tw' not in param_eq:
        if verbose:
            print("Wall temperature not defined. Using default Tw = 300 K.")
        param_eq['Tw'] = 300

    if len(env) > 1:
        param_eq = environment(param_eq, *env[:15])  # Atmospheric setup

    # Convert angles to radians
    aoa_rad = [radians(a) for a in aoa_deg]
    aos_rad = [radians(a) for a in aos_deg]

    # ----------------------------
    # Call main coefficient engine
    # ----------------------------
    output_path = calc_coeff(
        mod_path, res_path, aoa_rad, aos_rad,
        param_eq, flag_shadow, flag_solar, del_files,
        delete_temp_files=del_files, verbose=verbose
    )

    return output_path

