import numpy as np
from scipy.special import erf
from .ADBSatConstants import ConstantsData, EnvironmentData


def coeff_schaaf(param_eq, delta):
    """
    Computes aerodynamic coefficients using the Schaaf and Chambre model
    for a flat plate in free molecular flow (FMF).

    This semi-analytical model accounts for thermal accommodation and surface reflection.

    Parameters:
        param_eq (EnvironmentData or dict): Must contain:
            - sigmaN (float): Normal momentum accommodation coefficient [0–2]
            - sigmaT (float): Tangential momentum accommodation coefficient [0–1]
            - s (float): Speed ratio (freestream speed / thermal speed)
            - Tw (float): Wall temperature [K]
            - Tinf (float): Freestream temperature [K]
        delta (np.ndarray): Angles between surface normal and flow direction [rad]

    Returns:
        tuple of np.ndarray:
            - cp: pressure coefficient
            - ctau: shear stress coefficient
            - cd: drag coefficient
            - cl: lift coefficient
    """
    sigmaN = param_eq['sigmaN']
    sigmaT = param_eq['sigmaT']
    s = param_eq['s']
    Tw = param_eq['Tw']
    Tinf = param_eq['Tinf']

    cos_d = np.cos(delta)
    sin_d = np.sin(delta)

    # Pressure coefficient cp
    cp = (1 / s**2) * (
        ((2 - sigmaN) * s / np.sqrt(np.pi) * cos_d + sigmaN / 2 * np.sqrt(Tw / Tinf)) * 
        np.exp(-s**2 * cos_d**2) +
        ((2 - sigmaN) * (0.5 + s**2 * cos_d**2) +
         sigmaN / 2 * np.sqrt(Tw / Tinf) * np.sqrt(np.pi) * s * cos_d) *
        (1 + erf(s * cos_d))
    )

    # Shear stress coefficient ctau
    ctau = (sigmaT * sin_d / (s * np.sqrt(np.pi))) * (
        np.exp(-s**2 * cos_d**2) + s * np.sqrt(np.pi) * cos_d * (1 + erf(s * cos_d))
    )

    # Total drag and lift coefficients
    cd = cp * cos_d + ctau * sin_d
    cl = cp * sin_d - ctau * cos_d

    return cp, ctau, cd, cl

