import numpy as np
from .ADBSatConstants import ConstantsData, EnvironmentData


def coeff_newton(param_eq, delta):
    """
    Calculates aerodynamic coefficients using Newton's impact theory 
    for a flat plate in free molecular flow.

    This is a simplified high-speed model assuming inelastic particle-surface interaction
    with zero tangential shear forces.

    Parameters:
        param_eq (EnvironmentData or dict): Unused in this model, included for interface compatibility.
        delta (np.ndarray): Array of angles between surface normal and flow direction [rad].

    Returns:
        tuple of np.ndarray:
            - cp: pressure coefficient
            - ctau: shear stress coefficient (always zero)
            - cd: drag coefficient
            - cl: lift coefficient
    """
    # Pressure coefficient: cp = 2 * cos²(delta)
    cp = 2 * (np.cos(delta) ** 2)

    # No shear in Newtonian impact model
    ctau = np.zeros_like(delta)

    # Zero out coefficients for backward-facing surfaces
    cp[delta > np.pi / 2] = 0
    ctau[delta > np.pi / 2] = 0

    # Combine for total drag and lift
    cd = cp * np.cos(delta) + ctau * np.sin(delta)
    cl = cp * np.sin(delta) - ctau * np.cos(delta)

    return cp, ctau, cd, cl

