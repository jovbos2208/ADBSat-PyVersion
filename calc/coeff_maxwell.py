import numpy as np
from scipy.special import erf
from .ADBSatConstants import ConstantsData, EnvironmentData


def coeff_maxwell(param_eq, delta):
    """
    Computes aerodynamic coefficients using Maxwell's free molecular flow model.

    This model accounts for diffuse and specular reflection from the surface,
    incorporating thermal speed ratios and wall temperature effects.

    Parameters:
        param_eq (EnvironmentData or dict): Must contain:
            - alpha (float): Energy accommodation coefficient (0 = specular, 1 = diffuse)
            - Tw (float): Wall temperature [K]
            - Tinf (float): Freestream temperature [K]
            - s (float): Speed ratio (Vinf / most probable thermal speed)
        delta (np.ndarray): Array of angles between surface normal and flow direction [rad]

    Returns:
        tuple of np.ndarray:
            - cp: pressure coefficient
            - ctau: shear stress coefficient
            - cd: drag coefficient
            - cl: lift coefficient
    """
    const = ConstantsData()
    alpha = param_eq['alpha']
    Tw = param_eq['Tw']
    Tinf = param_eq['Tinf']
    s = param_eq['s']

    f = 1 - alpha
    theta = np.pi / 2 - delta  # angle between surface and flow (wall-aligned)

    # Pressure drag coefficient
    cd = 2 * ((1 - f * np.cos(2 * theta)) / (np.sqrt(np.pi) * s)) * np.exp(-s**2 * np.sin(theta)**2) + \
         (np.sin(theta) / s**2) * (1 + 2 * s**2 + f * (1 - 2 * s**2 * np.cos(2 * theta))) * erf(s * np.sin(theta)) + \
         ((1 - f) / s) * np.sqrt(np.pi) * np.sin(theta)**2 * np.sqrt(Tw / Tinf)

    # Lift coefficient
    cl = ((4 * f) / (np.sqrt(np.pi) * s)) * np.sin(theta) * np.cos(theta) * np.exp(-s**2 * np.sin(theta)**2) + \
         (np.cos(theta) / s**2) * (1 + f * (1 + 4 * s**2 * np.sin(theta)**2)) * erf(s * np.sin(theta)) + \
         ((1 - f) / s) * np.sqrt(np.pi) * np.sin(theta) * np.cos(theta) * np.sqrt(Tw / Tinf)

    # Set zero for back-facing panels
    cd[delta > np.pi / 2] = 0
    cl[delta > np.pi / 2] = 0

    # Convert to pressure and shear coefficients
    cp = cd * np.cos(delta) + cl * np.sin(delta)
    ctau = cd * np.sin(delta) - cl * np.cos(delta)

    return cp, ctau, cd, cl

