import numpy as np
from .ADBSatConstants import ConstantsData


def coeff_storchHyp(param_eq, delta):
    """
    Computes aerodynamic coefficients for a flat plate in a hyperthermal flow
    using Storch's model.

    The model accounts for both normal and tangential momentum accommodation and 
    uses a simplified kinetic formulation valid in the free molecular regime.

    Parameters:
        param_eq (dict): Dictionary with:
            - sigmaN (float): Normal momentum accommodation coefficient [0–2]
            - sigmaT (float): Tangential momentum accommodation coefficient [0–1]
            - Tw (float): Wall temperature [K]
            - Rmean (float): Mean specific gas constant [J/kg·K]
            - vinf (float): Freestream velocity [m/s]
        delta (np.ndarray): Array of angles between surface normal and freestream direction [rad]

    Returns:
        tuple of np.ndarray:
            - cp: pressure coefficient
            - ctau: shear stress coefficient
            - cd: drag coefficient
            - cl: lift coefficient
    """
    constants = ConstantsData()

    sigmaN = param_eq['sigmaN']
    sigmaT = param_eq['sigmaT']
    Tw = param_eq['Tw']
    Rmean = param_eq['Rmean']
    Vinf = param_eq['vinf']

    cos_d = np.cos(delta)
    sin_d = np.sin(delta)

    # Most probable thermal wall velocity
    Vw = np.sqrt((np.pi * Rmean * Tw) / 2)

    # Pressure coefficient
    cp = 2 * cos_d * (sigmaN * (Vw / Vinf) + (2 - sigmaN) * cos_d)

    # Shear stress coefficient
    ctau = 2 * cos_d * sin_d * sigmaT

    # Zero out coefficients for backward-facing surfaces
    cp[delta > np.pi / 2] = 0
    ctau[delta > np.pi / 2] = 0

    # Total drag and lift coefficients
    cd = cp * cos_d + ctau * sin_d
    cl = cp * sin_d - ctau * cos_d

    return cp, ctau, cd, cl

