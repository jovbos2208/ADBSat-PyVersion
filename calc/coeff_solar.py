import numpy as np


def coeff_solar(delta, param_eq):
    """
    Computes solar radiation pressure coefficients on a flat plate using 
    the Luthcke et al. (1997) model.

    Assumes no transmission (opaque surface) and:
        absorptivity + specular reflectivity + diffuse reflectivity = 1

    Parameters:
        delta (np.ndarray): Array of incident angles [rad] between solar vector
                            and surface normal.
        param_eq (dict): Dictionary with:
            - sol_cR (float): Specular reflectivity coefficient (rho)
            - sol_cD (float): Diffuse reflectivity coefficient (delta)

    Returns:
        tuple of np.ndarray:
            - cn: Normal coefficient (force component normal to the surface)
            - cs: Incident coefficient (force component in direction of sunlight)
    """
    rho = param_eq['sol_cR']  # Specular reflectivity
    d = param_eq['sol_cD']    # Diffuse reflectivity

    # Normal force component due to solar pressure
    cn = 2 * ((d / 3) * np.cos(delta) + rho * np.cos(delta) ** 2)

    # Incident (or projected) force component
    cs = (1 - rho) * np.cos(delta)
    cs[cs < 0] = 0  # No force from back-facing panels

    return cn, cs

