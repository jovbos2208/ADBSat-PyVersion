import numpy as np
from .ADBSatConstants import EnvironmentData


def coeff_cook(param_eq, delta):
    """
    Calculates aerodynamic coefficients for a flat plate using Cook's 
    hyperthermal free molecular flow model.

    Parameters:
        param_eq (EnvironmentData): Object or dict containing:
            - alpha (float): energy accommodation coefficient
            - Tw (float): wall temperature [K]
            - Rmean (float): mean specific gas constant [J/kg·K]
            - vinf (float): freestream velocity [m/s]
        delta (np.ndarray): Array of surface-normal-to-flow angles [radians].

    Returns:
        tuple of np.ndarray:
            - cp: pressure coefficient
            - ctau: shear stress coefficient
            - cd: drag coefficient
            - cl: lift coefficient
    """
    alpha = param_eq['alpha']
    Tw = param_eq['Tw']
    Rmean = param_eq['Rmean']
    vinf = param_eq['vinf']

    # Effective freestream temperature
    Tinf = vinf**2 / (3 * Rmean)

    # Compute drag and lift coefficients
    cd = 2.0 * np.cos(delta) * (
        1 + (2 / 3) * np.sqrt(1 + alpha * (Tw / Tinf - 1)) * np.cos(delta)
    )
    cl = (4 / 3) * np.sqrt(1 + alpha * (Tw / Tinf - 1)) * np.sin(delta) * np.cos(delta)

    # Set values to zero for back-facing panels (delta ≥ 90°)
    cd[delta >= np.pi / 2] = 0
    cl[delta >= np.pi / 2] = 0

    # Compute pressure and shear coefficients
    cp = cd * np.cos(delta) + cl * np.sin(delta)
    ctau = cd * np.sin(delta) - cl * np.cos(delta)

    return cp, ctau, cd, cl

