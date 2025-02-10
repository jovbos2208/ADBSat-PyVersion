import numpy as np

def coeff_solar(delta, param_eq):
    """
    Calculates solar coefficients for a flat plate using Luthcke et al 1997 formula.
    
    Adsorptivity (alpha) + Specular Reflectivity (rho) + Diffuse Reflectivity (delta) = 1
    Transmissivity = 0

    Parameters:
        delta (numpy.ndarray): Array of angles between the flow and the surface normal (radians).
        param_eq (dict): Dictionary containing:
            - sol_cR: Specular reflectivity component.
            - sol_cD: Diffuse reflectivity component.

    Returns:
        tuple: (cn, cs) where:
            - cn: Normal coefficient (numpy array).
            - cs: Incident coefficient (numpy array).
    """
    # Calculate normal coefficient
    cn = 2 * ((param_eq['sol_cD'] / 3) * np.cos(delta) + param_eq['sol_cR'] * np.cos(delta)**2)
    
    # Calculate incident coefficient
    cs = (1 - param_eq['sol_cR']) * np.cos(delta)
    cs[cs < 0] = 0  # Ensure cs is non-negative

    return cn, cs
