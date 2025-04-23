import numpy as np
from scipy.special import erf
from .ADBSatConstants import ConstantsData

def coeff_DRIA(param_eq, delta):
    const = ConstantsData()
    Tinf = param_eq['Tinf']
    Vinf = param_eq['vinf']
    alpha = param_eq['alpha']
    Tw = param_eq['Tw']
    gam = param_eq['gamma']  # typically cos(delta)
    ell = param_eq['ell']    # typically sin(delta)
    massConc = param_eq['massConc']

    # Molecular masses in g/mol
    molecular_masses = np.array([
        const.mHe, const.mO, const.mN2, const.mO2,
        const.mAr, const.mH, const.mN, const.mAnO
    ])

    n_species = len(molecular_masses)
    n_panels = len(delta)

    cp_j = np.zeros((n_species, n_panels))
    ctau_j = np.zeros((n_species, n_panels))

    for j, m_gmol in enumerate(molecular_masses):
        m_kg = m_gmol / 1000 / const.NA

        s = Vinf / np.sqrt(2 * const.kb * Tinf / m_kg)  # thermal speed ratio

        P = np.exp(-gam**2 * s**2) / s
        G = 1 / (2 * s**2)
        Q = 1 + G
        Z = 1 + erf(gam * s)

        R = const.R / m_gmol * 1e3  # specific gas constant in J/(kg K)
        Vratio = np.sqrt(0.5 * (1 + alpha * ((4 * R * Tw) / Vinf**2 - 1)))  # from Doornbos 2012

        cp_j[j, :] = P / np.sqrt(np.pi) + gam * Q * Z + 0.5 * gam * Vratio * (gam * np.sqrt(np.pi) * Z + P)
        ctau_j[j, :] = ell * G * Z + 0.5 * ell * Vratio * (gam * np.sqrt(np.pi) * Z + P)

    cp = np.sum(cp_j * massConc[:8, None], axis=0)
    ctau = np.sum(ctau_j * massConc[:8, None], axis=0)

    cd = cp * np.cos(delta) + ctau * np.sin(delta)
    cl = cp * np.sin(delta) - ctau * np.cos(delta)

    return cp, ctau, cd, cl
