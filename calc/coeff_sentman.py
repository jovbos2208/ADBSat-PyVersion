import numpy as np
from scipy.special import erf
from .ADBSatConstants import ConstantsData


def coeff_sentman(param_eq, delta):
    """
    Calculates aerodynamic coefficients using the Sentman model 
    for free molecular flow over a flat plate.

    This model accounts for thermal accommodation and wall temperature,
    and is based on kinetic theory approximations.

    Parameters:
        param_eq (dict): Dictionary containing:
            - alpha (float): Energy accommodation coefficient
            - s (float): Speed ratio (Vinf / thermal speed)
            - Tw (float): Wall temperature [K]
            - Tinf (float): Freestream temperature [K]
        delta (np.ndarray): Angles between surface normal and freestream direction [rad]

    Returns:
        tuple of np.ndarray:
            - cp: pressure coefficient
            - ctau: shear stress coefficient
            - cd: drag coefficient
            - cl: lift coefficient
    """
    constants = ConstantsData()

    alpha = param_eq['alpha']
    s = param_eq['s']
    Tw = param_eq['Tw']
    Tinf = param_eq['Tinf']

    # Intermediate wall temperature Ti (internal translation)
    Ti = 0.5 * s**2 * Tinf

    cos_d = np.cos(delta)
    sin_d = np.sin(delta)

    # Pressure coefficient cp
    cp = (
        (cos_d**2 + 1 / (2 * s**2)) * (1 + erf(s * cos_d)) +
        cos_d / (np.sqrt(np.pi) * s) * np.exp(-s**2 * cos_d**2) +
        0.5 * np.sqrt(0.5 * (1 + alpha * (Tw / Ti - 1))) *
        (np.sqrt(np.pi) * cos_d * (1 + erf(s * cos_d)) +
         1 / s * np.exp(-s**2 * cos_d**2))
    )

    # Shear stress coefficient ctau
    ctau = (
        sin_d * cos_d * (1 + erf(s * cos_d)) +
        sin_d / (s * np.sqrt(np.pi)) * np.exp(-s**2 * cos_d**2)
    )

    # Combine to total drag and lift
    cd = cp * cos_d + ctau * sin_d
    cl = cp * sin_d - ctau * cos_d

    return cp, ctau, cd, cl

