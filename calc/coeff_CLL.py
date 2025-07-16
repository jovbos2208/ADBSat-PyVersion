import numpy as np
from scipy.special import erf
from .ADBSatConstants import ConstantsData


def fitted_parameters(i, alphaN):
    """
    Computes the fitted CLL model parameters (beta, gamma, delta, zeta) 
    based on species index and the normal energy accommodation coefficient.

    Parameters:
        i (int): Index for the molecular species.
                 1 = He, 2 = O, 3 = N2, 4 = O2, 5 = H, 6 = N
        alphaN (numpy.ndarray): Normal energy accommodation coefficient
                                (vectorized across material surface).

    Returns:
        tuple of numpy.ndarray: Arrays of beta, gamma, delta, zeta values.
    """
    beta = np.ones_like(alphaN)
    gamma = np.ones_like(alphaN)
    delta = np.ones_like(alphaN)
    zeta = np.ones_like(alphaN)

    if i == 1:  # Helium
        beta[(alphaN > 0.95)] = 6.2
        gamma[(alphaN > 0.95)] = 0.38
        delta[(alphaN > 0.95)] = 3.3
        zeta[(alphaN > 0.95)] = 0.74

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

    elif i == 2:  # Oxygen (O)
        beta *= 5.85
        gamma *= 0.2
        delta *= 0.48
        zeta *= 31.0

    elif i == 3:  # Nitrogen (N2)
        beta *= 6.6
        gamma *= 0.22
        delta *= 0.48
        zeta *= 35.0

    elif i == 4:  # Oxygen molecule (O2)
        beta *= 6.3
        gamma *= 0.26
        delta *= 0.42
        zeta *= 20.5

    elif i == 5:  # Hydrogen (H)
        beta[(alphaN > 0.95)] = 3.9
        gamma[(alphaN > 0.95)] = 0.195
        delta[(alphaN > 0.95)] = 1.4
        zeta[(alphaN > 0.95)] = 0.3

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

    elif i == 6:  # Atomic nitrogen (N)
        beta *= 4.9
        gamma *= 0.32
        delta *= 0.42
        zeta *= 8.0

    return beta, gamma, delta, zeta


def coeff_cll(env_data, delta):
    """
    Computes aerodynamic coefficients using the Cercignani-Lampis-Lord (CLL) model.

    Parameters:
        env_data (dict): Environmental data with keys:
            - alphaN (array): normal accommodation coefficients
            - sigmaT (array): tangential momentum accommodation
            - vinf (float): freestream velocity
            - Tw (float): wall temperature
            - Tinf (float): freestream temperature
            - rho (array): species densities (at least indices 0-9 expected)
        delta (array): angles between surface normals and freestream direction [rad]

    Returns:
        tuple of arrays: cp, ctau, cd, cl
            - cp: pressure coefficient
            - ctau: shear stress coefficient
            - cd: drag coefficient
            - cl: lift coefficient
    """
    constants = ConstantsData()

    # Molecular masses for He, O, N2, O2, H, N [g/mol]
    M_j = np.array([constants.mHe, constants.mO, constants.mN2,
                    constants.mO2, constants.mH, constants.mN])

    alphaN = env_data['alphaN']
    sigmaT = env_data['sigmaT']
    vinf = env_data['vinf']
    Tw = env_data['Tw']
    Tinf = env_data['Tinf']
    rho = np.copy(env_data['rho'])

    # Remove species not modeled (e.g., Ar, AO)
    rho[4] = 0
    rho[7] = 0
    rho[8] = 0
    rho[9] = 0
    rho = rho[rho != 0]

    Ns = len(rho)
    Nelem = len(delta)

    cp_j = np.zeros((Ns, Nelem))
    ctau_j = np.zeros((Ns, Nelem))
    cd_j = np.zeros((Ns, Nelem))
    cl_j = np.zeros((Ns, Nelem))

    for j in range(Nelem):
        for i in range(Ns):
            s_i = vinf / np.sqrt(2 * (constants.kb / ((M_j[i] / constants.NA) / 1000) * Tinf))

            beta_fp, gamma_fp, delta_fp, zeta_fp = fitted_parameters(i + 1, alphaN)

            x_var = s_i * np.cos(delta[j])
            gamma1 = (1 / np.sqrt(np.pi)) * (
                x_var * np.exp(-x_var**2) + 
                (np.sqrt(np.pi) / 2) * (1 + 2 * x_var**2) * (1 + erf(x_var))
            )
            gamma2 = (1 / np.sqrt(np.pi)) * (
                np.exp(-x_var**2) + np.sqrt(np.pi) * x_var * (1 + erf(x_var))
            )

            if alphaN[j] < 1:
                cp_j[i, j] = (1 / s_i**2) * (
                    (1 + np.sqrt(1 - alphaN[j])) * gamma1 + 
                    0.5 * np.exp(-beta_fp[j] * (1 - alphaN[j])**gamma_fp[j]) *
                    (Tw / Tinf)**delta_fp[j] * (zeta_fp[j] / s_i) *
                    np.sqrt(Tw / Tinf) * np.sqrt(np.pi) * gamma2
                )
                ctau_j[i, j] = (sigmaT[j] * np.sin(delta[j])) / s_i * gamma2
            else:
                cp_j[i, j] = (1 / s_i**2) * (
                    ((2 - alphaN[j]) * s_i / np.sqrt(np.pi) * np.cos(delta[j]) +
                     alphaN[j] / 2 * (Tw / Tinf)**0.5) * np.exp(-s_i**2 * np.cos(delta[j])**2) +
                    ((2 - alphaN[j]) * (0.5 + s_i**2 * np.cos(delta[j])**2) +
                     alphaN[j] / 2 * (Tw / Tinf)**0.5 * np.sqrt(np.pi) * s_i * np.cos(delta[j])) *
                    (1 + erf(s_i * np.cos(delta[j])))
                )
                ctau_j[i, j] = (sigmaT[j] * np.sin(delta[j])) / (s_i * np.sqrt(np.pi)) * (
                    np.exp(-s_i**2 * np.cos(delta[j])**2) +
                    s_i * np.sqrt(np.pi) * np.cos(delta[j]) * (1 + erf(s_i * np.cos(delta[j])))
                )

            cd_j[i, j] = cp_j[i, j] * np.cos(delta[j]) + ctau_j[i, j] * np.sin(delta[j])
            cl_j[i, j] = cp_j[i, j] * np.sin(delta[j]) - ctau_j[i, j] * np.cos(delta[j])

    xi_j = rho / np.sum(rho)
    m_avg = np.sum(xi_j * M_j[:Ns] / constants.NA * rho)

    sum_cp = np.sum(xi_j[:, None] * M_j[:Ns, None] / constants.NA * rho[:, None] * cp_j, axis=0)
    sum_ctau = np.sum(xi_j[:, None] * M_j[:Ns, None] / constants.NA * rho[:, None] * ctau_j, axis=0)
    sum_cd = np.sum(xi_j[:, None] * M_j[:Ns, None] / constants.NA * rho[:, None] * cd_j, axis=0)
    sum_cl = np.sum(xi_j[:, None] * M_j[:Ns, None] / constants.NA * rho[:, None] * cl_j, axis=0)

    cp = sum_cp / m_avg
    ctau = sum_ctau / m_avg
    cd = sum_cd / m_avg
    cl = sum_cl / m_avg

    return cp, ctau, cd, cl

