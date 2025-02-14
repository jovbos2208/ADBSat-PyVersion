import numpy as np
from scipy.special import erf
from .ADBSatConstants import ConstantsData
def fitted_parameters(i, alphaN):
    """
    Computes fitted parameters for the CLL model according to Walker et al. (2014).
    
    Parameters:
    i : int
        Index representing the molecular species (1 = He, 2 = O, etc.).
    alphaN : numpy array
        Normal accommodation coefficient (vectorized for multiple materials).

    Returns:
    beta, gamma, delta, zeta : numpy arrays
        Fitted parameters for the CLL model.
    """

    # Initialisiere Arrays mit Einsen (wie `temp = ones(1,length(alphaN))` in MATLAB)
    beta = np.ones_like(alphaN)
    gamma = np.ones_like(alphaN)
    delta = np.ones_like(alphaN)
    zeta = np.ones_like(alphaN)

    if i == 1:  # Helium (He)
        beta[(alphaN > 0.95) & (alphaN <= 1)] = 6.2
        gamma[(alphaN > 0.95) & (alphaN <= 1)] = 0.38
        delta[(alphaN > 0.95) & (alphaN <= 1)] = 3.3
        zeta[(alphaN > 0.95) & (alphaN <= 1)] = 0.74

        beta[(alphaN > 0.90) & (alphaN <= 0.95)] = 3.8
        gamma[(alphaN > 0.90) & (alphaN <= 0.95)] = 0.52
        delta[(alphaN > 0.90) & (alphaN <= 0.95)] = 3.4
        zeta[(alphaN > 0.90) & (alphaN <= 0.95)] = 1.12

        beta[(alphaN > 0.50) & (alphaN <= 0.90)] = 3.45
        gamma[(alphaN > 0.50) & (alphaN <= 0.90)] = 0.52
        delta[(alphaN > 0.50) & (alphaN <= 0.90)] = 2.4
        zeta[(alphaN > 0.50) & (alphaN <= 0.90)] = 0.93

        beta[(alphaN > 0) & (alphaN <= 0.50)] = 0.08
        gamma[(alphaN > 0) & (alphaN <= 0.50)] = 0.52
        delta[(alphaN > 0) & (alphaN <= 0.50)] = 4.2
        zeta[(alphaN > 0) & (alphaN <= 0.50)] = 1.1

    elif i == 2:  # Sauerstoff (O)
        beta *= 5.85
        gamma *= 0.2
        delta *= 0.48
        zeta *= 31.0

    elif i == 3:  # Stickstoff (N2)
        beta *= 6.6
        gamma *= 0.22
        delta *= 0.48
        zeta *= 35.0

    elif i == 4:  # Sauerstoff (O2)
        beta *= 6.3
        gamma *= 0.26
        delta *= 0.42
        zeta *= 20.5

    elif i == 5:  # Wasserstoff (H)
        beta[(alphaN > 0.95) & (alphaN <= 1)] = 3.9
        gamma[(alphaN > 0.95) & (alphaN <= 1)] = 0.195
        delta[(alphaN > 0.95) & (alphaN <= 1)] = 1.4
        zeta[(alphaN > 0.95) & (alphaN <= 1)] = 0.3

        beta[(alphaN > 0.90) & (alphaN <= 0.95)] = 3.5
        gamma[(alphaN > 0.90) & (alphaN <= 0.95)] = 0.42
        delta[(alphaN > 0.90) & (alphaN <= 0.95)] = 2.0
        zeta[(alphaN > 0.90) & (alphaN <= 0.95)] = 0.72

        beta[(alphaN > 0.50) & (alphaN <= 0.90)] = 3.45
        gamma[(alphaN > 0.50) & (alphaN <= 0.90)] = 0.52
        delta[(alphaN > 0.50) & (alphaN <= 0.90)] = 2.4
        zeta[(alphaN > 0.50) & (alphaN <= 0.90)] = 0.93

        beta[(alphaN > 0) & (alphaN <= 0.50)] = 0.095
        gamma[(alphaN > 0) & (alphaN <= 0.50)] = 0.465
        delta[(alphaN > 0) & (alphaN <= 0.50)] = 2.9
        zeta[(alphaN > 0) & (alphaN <= 0.50)] = 0.92

    elif i == 6:  # Stickstoff (N)
        beta *= 4.9
        gamma *= 0.32
        delta *= 0.42
        zeta *= 8.0

    return beta, gamma, delta, zeta



import numpy as np
from scipy.special import erf


def coeff_cll(env_data, delta):
    """
    Computes aerodynamic coefficients using the CLL model.

    Parameters:
    env_data : dict
        Contains parameters such as alphaN, sigmaT, Tw, Tinf, vinf, and species densities.
    delta : numpy array
        Angles between the surface normal and the flow.

    Returns:
    cp : numpy array of shape (Nelem,)
        Normal pressure coefficient.
    ctau : numpy array of shape (Nelem,)
        Shear stress coefficient.
    cd : numpy array of shape (Nelem,)
        Drag coefficient.
    cl : numpy array of shape (Nelem,)
        Lift coefficient.
    """

    constants = ConstantsData()
    
    # Molekülmassen [g/mol] für {He, O, N2, O2, H, N}
    M_j = np.array([
        constants.mHe, constants.mO, constants.mN2,
        constants.mO2, constants.mH, constants.mN
    ])

    # Parameter aus `env_data`
    alphaN = env_data['alphaN']  # (array der Länge Ns)
    sigmaT = env_data['sigmaT']
    vinf = env_data['vinf']
    Tw = env_data['Tw']
    Tinf = env_data['Tinf']
    rho = np.copy(env_data['rho'])  # Speziesdichten (Kopie erstellen, um Original nicht zu verändern)

    # Entferne Nullwerte aus rho, Ignoriere Ar,Total und AO
    rho[4] = 0 
    rho[7] = 0
    rho[8] = 0
    rho[9] = 0
    rho = rho[rho != 0]

    # Anzahl der Spezies und Anzahl der Elemente
    Ns = len(rho)
    Nelem = len(delta)

    # Initialisierung der Gesamtkraft-Koeffizienten als Arrays der Länge `Nelem`
    ctau_j = np.zeros((Ns, Nelem))
    cp_j = np.zeros((Ns, Nelem))
    cd_j = np.zeros((Ns, Nelem))
    cl_j = np.zeros((Ns, Nelem))

    # Berechnung für jedes Element in `delta`
    for j in range(Nelem):
        for i in range(Ns):
            s_i = vinf / np.sqrt(2 * (constants.kb / ((M_j[i] / constants.NA) / 1000) * Tinf))

            # Berechne angepasste Parameter nach Walker et al.
            beta_fp, gamma_fp, delta_fp, zeta_fp = fitted_parameters(i, alphaN)

            # Hilfsfunktionen für Gamma-Werte
            x_var = s_i * np.cos(delta[j])
            gamma1 = (1 / np.sqrt(np.pi)) * (x_var * np.exp(-x_var**2) + np.sqrt(np.pi) / 2 * (1 + 2 * x_var**2) * (1 + erf(x_var)))
            gamma2 = (1 / np.sqrt(np.pi)) * (np.exp(-x_var**2) + np.sqrt(np.pi) * x_var * (1 + erf(x_var)))

            # Berechnung der aerodynamischen Koeffizienten
            if alphaN[j] < 1:
                cp_j[i, j] = (1 / s_i**2) * ((1 + np.sqrt(1 - alphaN[j])) * gamma1 + 
                            0.5 * (np.exp(-beta_fp * (1 - alphaN[j])**gamma_fp) * (Tw / Tinf)**delta_fp * (zeta_fp / s_i)) *
                            np.sqrt(Tw / Tinf) * np.sqrt(np.pi) * gamma2)
                ctau_j[i, j] = (sigmaT[j] * np.sin(delta[j])) / s_i * gamma2
            else:
                cp_j[i, j] = (1 / s_i**2) * (((2 - alphaN[j]) * s_i / np.sqrt(np.pi) * np.cos(delta[j]) + alphaN[j] / 2 * (Tw / Tinf)**0.5) * np.exp(-s_i**2 * np.cos(delta[j])**2) +
                            ((2 - alphaN[j]) * (0.5 + s_i**2 * np.cos(delta[j])**2) + alphaN[j] / 2 * (Tw / Tinf)**0.5 * np.sqrt(np.pi) * s_i * np.cos(delta[j])) * (1 + erf(s_i * np.cos(delta[j]))))
                ctau_j[i, j] = (sigmaT[j] * np.sin(delta[j])) / (s_i * np.sqrt(np.pi)) * (np.exp(-s_i**2 * np.cos(delta[j])**2) + s_i * np.sqrt(np.pi) * np.cos(delta[j]) * (1 + erf(s_i * np.cos(delta[j]))))

            # Berechnung von cd und cl für diese Spezies
            cd_j[i, j] = cp_j[i, j] * np.cos(delta[j]) + ctau_j[i, j] * np.sin(delta[j])
            cl_j[i, j] = cp_j[i, j] * np.sin(delta[j]) - ctau_j[i, j] * np.cos(delta[j])

    # Berechnung der Gesamtkoeffizienten für jedes Panel
    xi_j = rho / np.sum(rho)  # Molenbruch der Spezies
    m_avg = np.sum(xi_j * M_j[:Ns] / constants.NA * rho)  # Durchschnittliche Masse des Gasgemischs

    sum_cp = np.sum(xi_j[:, np.newaxis] * M_j[:Ns, np.newaxis] / constants.NA * rho[:, np.newaxis] * cp_j, axis=0)
    sum_ctau = np.sum(xi_j[:, np.newaxis] * M_j[:Ns, np.newaxis] / constants.NA * rho[:, np.newaxis] * ctau_j, axis=0)
    sum_cd = np.sum(xi_j[:, np.newaxis] * M_j[:Ns, np.newaxis] / constants.NA * rho[:, np.newaxis] * cd_j, axis=0)
    sum_cl = np.sum(xi_j[:, np.newaxis] * M_j[:Ns, np.newaxis] / constants.NA * rho[:, np.newaxis] * cl_j, axis=0)

    # Normierung durch `m_avg`
    cp = sum_cp / m_avg
    ctau = sum_ctau / m_avg
    cd = sum_cd / m_avg
    cl = sum_cl / m_avg

    return cp, ctau, cd, cl
