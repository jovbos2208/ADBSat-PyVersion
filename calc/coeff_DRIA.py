import numpy as np
from scipy.special import erf
from .ADBSatConstants import ConstantsData


def coeff_DRIA(param_eq, delta):
    """
    Computes aerodynamic coefficients using the DRIA (Doornbos Refined Impact Approximation) model.

    This model accounts for molecular species with mass-weighted contributions and surface interaction.

    Parameters:
        param_eq (dict): Dictionary containing:
            - Tinf (float): Freestream temperature [K]
            - vinf (float): Freestream velocity [m/s]
            - alpha (float): Energy accommodation coefficient
            - Tw (float): Wall temperature [K]
            - gamma (np.ndarray): Cosine of angle delta (precomputed)
            - ell (np.ndarray): Sine of angle delta (precomputed)
            - massConc (np.ndarray): Species mass concentrations (shape: [8,])

        delta (np.ndarray): Angles between surface normal and flow direction [radians]

    Returns:
        tuple of np.ndarray:
            - cp: pressure coefficient
            - ctau: shear stress coefficient
            - cd: drag coefficient
            - cl: lift coefficient
    """
    const = ConstantsData()
    Tinf = param_eq['Tinf']
    Vinf = param_eq['vinf']
    alpha = param_eq['alpha']
    Tw = param_eq['Tw']
    gam = param_eq['gamma']
    ell = param_eq['ell']
    massConc = param_eq['massConc']

    # Molecular masses [g/mol] for 8 species
    molecular_masses = np.array([
        const.mHe, const.mO, const.mN2, const.mO2,
        const.mAr, const.mH, const.mN, const.mAnO
    ])

    n_species = len(molecular_masses)
    n_panels = len(delta)

    cp_j = np.zeros((n_species, n_panels))
    ctau_j = np.zeros((n_species, n_panels))

    for j, m_gmol in enumerate(molecular_masses):
        m_kg = m_gmol / 1000 / const.NA  # Convert to kg

        # Speed ratio
        s = Vinf / np.sqrt(2 * const.kb * Tinf / m_kg)

        # DRIA intermediate terms
        P = np.exp(-gam**2 * s**2) / s
        G = 1 / (2 * s**2)
        Q = 1 + G
        Z = 1 + erf(gam * s)

        # Species-specific gas constant [J/kg·K]
        R = const.R / m_gmol * 1e3

        # Velocity ratio term from Doornbos 2012
        Vratio = np.sqrt(0.5 * (1 + alpha * ((4 * R * Tw) / Vinf**2 - 1)))

        # Species contribution to cp and ctau
        cp_j[j, :] = P / np.sqrt(np.pi) + gam * Q * Z + 0.5 * gam * Vratio * (
            gam * np.sqrt(np.pi) * Z + P
        )
        ctau_j[j, :] = ell * G * Z + 0.5 * ell * Vratio * (
            gam * np.sqrt(np.pi) * Z + P
        )

    # Mass-weighted summation across all species
    cp = np.sum(cp_j * massConc[:8, None], axis=0)
    ctau = np.sum(ctau_j * massConc[:8, None], axis=0)

    # Drag and lift coefficients
    cd = cp * np.cos(delta) + ctau * np.sin(delta)
    cl = cp * np.sin(delta) - ctau * np.cos(delta)

    return cp, ctau, cd, cl

