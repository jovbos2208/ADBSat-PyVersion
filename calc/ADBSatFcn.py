import os
import numpy as np
from math import radians
from .ADBSatConstants import EnvironmentData, ConstantsData
from .environment import environment
from .calc_coeff import calc_coeff

def adbsat_fcn(mod_path, res_path, param_eq, aoa_deg, aos_deg, flag_shadow, flag_solar, env, del_files, verbose):
    """
    Main function for creating aerodynamic and solar coefficients files.

    Parameters:
        mod_path (str): Path of the input model file.
        res_path (str): Path for the output file.
        param_eq (dict): Parameters for GSI and SRP models.
        aoa_deg (list): List of angles of attack in degrees.
        aos_deg (list): List of angles of sideslip in degrees.
        flag_shadow (bool): Perform shadow analysis.
        flag_solar (bool): Calculate solar coefficients.
        env (list): Environmental parameters.
        del_files (bool): Flag for individual file cleanup.
        verbose (bool): Flag for verbose output.

    Returns:
        str: Path to the generated output file(s).
    """
    # Validate GSI model inputs
    gsi_model = param_eq.get('gsi_model', '').lower()
    if gsi_model == 'cll':
        if not all(key in param_eq for key in ['alphaN', 'sigmaT']):
            raise ValueError("CLL model requires 'alphaN' and 'sigmaT' inputs.")
    elif gsi_model in ['schaaf', 'storchhyp']:
        if not all(key in param_eq for key in ['sigmaN', 'sigmaT']):
            raise ValueError(f"{gsi_model} model requires 'sigmaN' and 'sigmaT' inputs.")
    elif gsi_model in ['cook', 'sentman', 'maxwell', 'dria']:
        if 'alpha' not in param_eq:
            raise ValueError(f"{gsi_model} model requires 'alpha' input.")
    elif gsi_model == 'newton':
        pass  # No specific parameters required
    else:
        raise ValueError("Unrecognized GSI model input string.")

    # Validate solar reflectivity inputs
    if flag_solar and not all(key in param_eq for key in ['sol_cR', 'sol_cD']):
        raise ValueError("SRP interaction model requires 'sol_cR' and 'sol_cD' inputs.")

    # Default wall temperature if not defined
    if 'Tw' not in param_eq:
        if verbose:
            print("Wall temperature not defined. Using default value of 300K.")
        param_eq['Tw'] = 300

    # Perform environmental calculations if necessary
    if len(env) > 1:
        param_eq = environment(param_eq, *env[:15])

    # Convert angles to radians
    aoa_rad = [radians(angle) for angle in aoa_deg]
    aos_rad = [radians(angle) for angle in aos_deg]

    # Perform coefficient calculations
    output_path = calc_coeff(mod_path, res_path, aoa_rad, aos_rad, param_eq, flag_shadow, flag_solar, del_files, verbose)

    return output_path
