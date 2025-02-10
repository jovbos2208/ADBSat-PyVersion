import numpy as np
from .ADBSatConstants import ConstantsData, EnvironmentData

def coeff_newton(param_eq, delta):
    """
    Calculates aerodynamic coefficients for a flat plate using Newton's formula.
    
    Parameters:
        param_eq (EnvironmentData): Structure containing necessary parameters.
        delta (ndarray): Array of angles between the surface normal and the flow (radians).
    
    Returns:
        tuple: cp, ctau, cd, cl (numpy arrays of coefficients for each panel).
    """
    cp = 2 * (np.cos(delta) ** 2)
    ctau = np.zeros_like(delta)

    # Backward-facing panels
    cp[delta > np.pi / 2] = 0
    ctau[delta > np.pi / 2] = 0

    cd = cp * np.cos(delta) + ctau * np.sin(delta)
    cl = cp * np.sin(delta) - ctau * np.cos(delta)

    return cp, ctau, cd, cl
