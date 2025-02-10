import numpy as np
from scipy.special import erf
from .ADBSatConstants import ConstantsData, EnvironmentData

def coeff_schaaf(param_eq, delta):
    """
    Calculates aerodynamic coefficients for a flat plate using Schaaf and Chambre formula.
    
    Parameters:
        param_eq (EnvironmentData): Structure containing necessary parameters.
        delta (ndarray): Array of angles between the surface normal and the flow (radians).
    
    Returns:
        tuple: cp, ctau, cd, cl (numpy arrays of coefficients for each panel).
    """
    sigmaN = param_eq['sigmaN']
    sigmaT = param_eq['sigmaT']
    s = param_eq['s']
    Tw = param_eq['Tw']
    Tinf = param_eq['Tinf']

    cp = (1 / s**2) * (
        ((2 - sigmaN) * s / np.sqrt(np.pi) * np.cos(delta) + sigmaN / 2 * (Tw / Tinf)**0.5) 
        * np.exp(-s**2 * (np.cos(delta))**2) +
        ((2 - sigmaN) * (0.5 + s**2 * (np.cos(delta))**2) + sigmaN / 2 * (Tw / Tinf)**0.5 * np.sqrt(np.pi) * s * np.cos(delta)) 
        * (1 + erf(s * np.cos(delta)))
    )

    ctau = (sigmaT * np.sin(delta) / (s * np.sqrt(np.pi))) * (
        np.exp(-s**2 * (np.cos(delta))**2) + s * np.sqrt(np.pi) * np.cos(delta) * (1 + erf(s * np.cos(delta)))
    )

    cd = cp * np.cos(delta) + ctau * np.sin(delta)
    cl = cp * np.sin(delta) - ctau * np.cos(delta)

    return cp, ctau, cd, cl
