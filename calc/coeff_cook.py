import numpy as np
from .ADBSatConstants import EnvironmentData

def coeff_cook(param_eq, delta):
    """
    Calculates "Hyperthermal" FMF coefficients for a flat plate using Cook's formula.

    Parameters:
        param_eq (EnvironmentData): Object containing environmental parameters.
        delta (numpy.ndarray): Array of angles between the surface normal and the flow [radians].

    Returns:
        tuple: cp, ctau, cd, cl (numpy arrays of coefficients for each panel).
    """
    alpha = param_eq['alpha']
    Tw = param_eq['Tw']
    Rmean = param_eq['Rmean']
    vinf = param_eq['vinf']
    
    # Compute kinetic temperature
    Tinf = vinf**2 / (3 * Rmean)

    # Calculate drag and lift coefficients
    cd = 2.0 * np.cos(delta) * (1 + (2/3) * np.sqrt(1 + alpha * (Tw / Tinf - 1)) * np.cos(delta))
    cl = (4/3) * np.sqrt(1 + alpha * (Tw / Tinf - 1)) * np.sin(delta) * np.cos(delta)

    # Set coefficients to zero for angles >= pi/2
    cd[delta >= np.pi / 2] = 0
    cl[delta >= np.pi / 2] = 0

    # Compute pressure and shear stress coefficients
    cp = cd * np.cos(delta) + cl * np.sin(delta)
    ctau = cd * np.sin(delta) - cl * np.cos(delta)

    return cp, ctau, cd, cl
