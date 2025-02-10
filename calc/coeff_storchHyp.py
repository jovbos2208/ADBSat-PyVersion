import numpy as np
from .ADBSatConstants import ConstantsData

def coeff_storchHyp(param_eq, delta):
    """
    Calculates hyperthermal aerodynamic coefficients for a flat plate using Storch's formula.

    Parameters:
        param_eq (dict): Dictionary containing model parameters (sigmaN, sigmaT, Tw, Rmean, Vinf).
        delta (numpy.ndarray): Array of angles between the surface normal and the flow (radians).

    Returns:
        tuple: cp, ctau, cd, cl (numpy arrays of coefficients for each panel).
    """
    constants = ConstantsData()

    sigmaN = param_eq['sigmaN']
    sigmaT = param_eq['sigmaT']
    Tw = param_eq['Tw']
    Rmean = param_eq['Rmean']
    Vinf = param_eq['vinf']

    Vw = np.sqrt((np.pi * Rmean * Tw) / 2)

    cp = 2 * np.cos(delta) * (sigmaN * (Vw / Vinf) + (2 - sigmaN) * np.cos(delta))
    ctau = 2 * np.cos(delta) * np.sin(delta) * sigmaT

    # Ensure coefficients are zero for backward-facing panels
    cp[delta > np.pi / 2] = 0
    ctau[delta > np.pi / 2] = 0

    cd = cp * np.cos(delta) + ctau * np.sin(delta)
    cl = cp * np.sin(delta) - ctau * np.cos(delta)

    return cp, ctau, cd, cl
