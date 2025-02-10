import numpy as np
from scipy.special import erf
from .ADBSatConstants import ConstantsData

def coeff_sentman(param_eq, delta):
    """
    Calculates aerodynamic coefficients for a flat plate using Sentman's formula.

    Parameters:
        param_eq (dict): Dictionary containing model parameters (alpha, s, Tw, Tatm).
        delta (numpy.ndarray): Array of angles between the surface normal and the flow (radians).

    Returns:
        tuple: cp, ctau, cd, cl (numpy arrays of coefficients for each panel).
    """
    constants = ConstantsData()

    alpha = param_eq['alpha']
    s = param_eq['s']
    Tw = param_eq['Tw']
    Tinf = param_eq['Tinf']

    Ti = 0.5 * s**2 * Tinf

    cp = ((np.cos(delta)**2 + 1 / (2 * s**2)) * (1 + erf(s * np.cos(delta))) +
          np.cos(delta) / (np.sqrt(np.pi) * s) * np.exp(-s**2 * np.cos(delta)**2) +
          0.5 * np.sqrt(0.5 * (1 + alpha * (Tw / Ti - 1))) *
          (np.sqrt(np.pi) * np.cos(delta) * (1 + erf(s * np.cos(delta))) +
           1 / s * np.exp(-s**2 * np.cos(delta)**2)))

    ctau = (np.sin(delta) * np.cos(delta) * (1 + erf(s * np.cos(delta))) +
            np.sin(delta) / (s * np.sqrt(np.pi)) * np.exp(-s**2 * np.cos(delta)**2))

    cd = cp * np.cos(delta) + ctau * np.sin(delta)
    cl = cp * np.sin(delta) - ctau * np.cos(delta)

    return cp, ctau, cd, cl
