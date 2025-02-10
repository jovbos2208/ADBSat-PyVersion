import numpy as np
from scipy.special import erf
from .ADBSatConstants import ConstantsData, EnvironmentData

def coeff_maxwell(param_eq, delta):
    """
    Calculates FMF aerodynamic coefficients according to Maxwell's model.
    
    Parameters:
        param_eq (EnvironmentData): Structure containing necessary parameters.
        delta (ndarray): Array of angles between the surface normal and the flow (radians).
    
    Returns:
        tuple: cp, ctau, cd, cl (numpy arrays of coefficients for each panel).
    """
    const = ConstantsData()
    alpha = param_eq['alpha']
    Tw = param_eq['Tw']
    Tinf = param_eq['Tinf']
    s = param_eq['s']

    f = 1 - alpha
    theta = np.pi / 2 - delta

    cd = 2 * ((1 - f * np.cos(2 * theta)) / (np.sqrt(np.pi) * s)) * np.exp(-s**2 * np.sin(theta)**2) + \
         (np.sin(theta) / s**2) * (1 + 2 * s**2 + f * (1 - 2 * s**2 * np.cos(2 * theta))) * \
         erf(s * np.sin(theta)) + \
         ((1 - f) / s) * np.sqrt(np.pi) * np.sin(theta)**2 * np.sqrt(Tw / Tinf)

    cl = ((4 * f) / (np.sqrt(np.pi) * s)) * np.sin(theta) * np.cos(theta) * np.exp(-s**2 * np.sin(theta)**2) + \
         (np.cos(theta) / s**2) * (1 + f * (1 + 4 * s**2 * np.sin(theta)**2)) * erf(s * np.sin(theta)) + \
         ((1 - f) / s) * np.sqrt(np.pi) * np.sin(theta) * np.cos(theta) * np.sqrt(Tw / Tinf)

    # Zero coefficients for panels facing backward
    cd[delta > np.pi / 2] = 0
    cl[delta > np.pi / 2] = 0

    cp = cd * np.cos(delta) + cl * np.sin(delta)
    ctau = cd * np.sin(delta) - cl * np.cos(delta)

    return cp, ctau, cd, cl
