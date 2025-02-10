import numpy as np
from scipy.special import erf
from .ADBSatConstants import ConstantsData, EnvironmentData

def coeff_DRIA(param_eq, delta):
    """
    Calculates aerodynamic coefficients for a flat plate using DRIA (Diffuse Reflection with Incomplete Accommodation).
    
    Parameters:
        param_eq (EnvironmentData): Structure containing necessary parameters.
        delta (ndarray): Array of angles between the surface normal and the flow (radians).
    
    Returns:
        tuple: cp, ctau, cd, cl (numpy arrays of coefficients for each panel).
    """
    const = ConstantsData()
    Tinf = param_eq['Tinf']
    Vinf = param_eq['vinf']
    alpha = param_eq['alpha']
    Tw = param_eq['Tw']
    massConc = param_eq['massConc']  

    # Molecular masses of species
    molecular_masses = np.array([const.mHe, const.mO, const.mN2, const.mO2, const.mAr, const.mH, const.mN, const.mAnO])

    gam = param_eq['gamma']
    ell = param_eq['ell']

    n_panels = len(delta)
    n_species = molecular_masses.shape[0]
    cp_j = np.zeros((n_species, n_panels))
    ctau_j = np.zeros((n_species, n_panels))

    for j, mass in enumerate(molecular_masses):
        s = Vinf / np.sqrt(2 * (const.kb / (mass / const.NA / 1000) * Tinf))
        

        P = np.exp(-gam**2 * s**2) / s
        G = 1 / (2 * s**2)
        Q = 1 + G
        Z = 1 + erf(gam * s)
        R = const.R / mass * 1e3  # Specific gas constant
        Vratio = np.sqrt(0.5 * (1 + alpha * ((4 * R * Tw) / Vinf**2 - 1)))

        cp_j[j, :] = P / np.sqrt(np.pi) + gam * Q * Z + 0.5 * gam * Vratio * (gam * np.sqrt(np.pi) * Z + P)
        ctau_j[j, :] = ell * G * Z + 0.5 * ell * Vratio * (gam * np.sqrt(np.pi) * Z + P)

    # Average coefficients
    cp = np.sum(cp_j * massConc[:8, None], axis=0)
    ctau = np.sum(ctau_j * massConc[:8, None], axis=0)

    cd = cp * np.cos(delta) + ctau * np.sin(delta)
    cl = cp * np.sin(delta) - ctau * np.cos(delta)

    return cp, ctau, cd, cl
