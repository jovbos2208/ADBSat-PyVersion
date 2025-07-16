import numpy as np
from .ADBSatConstants import ConstantsData, EnvironmentData
import pandas as pd


def environment(param_eq, database, idx, h):
    """
    Enriches the GSI model input parameters with atmospheric data using MSISE-00 model output.

    Assumes a circular orbit without co-rotation or wind. Uses data from a precomputed MSISE dataset.

    Parameters:
        param_eq (dict): Parameter dictionary to update (e.g. with vinf, Rmean, Tinf, etc.).
        database (pd.DataFrame): DataFrame of atmospheric values (MSISE-00 output).
        idx (int): Index into the database for the current simulation step.
        h (float): Altitude [m].

    Returns:
        dict: Updated `param_eq` dictionary with computed atmospheric values.
    """
    constants = ConstantsData()

    # Convert altitude to kilometers for consistency (e.g. with MSISE)
    alt_km = h / 1000.0

    # Select atmospheric values at given index
    atmosphere = database.iloc[idx].to_numpy()

    # Extract species densities into rho vector
    rho = np.array([
        atmosphere[4],  # He
        atmosphere[3],  # O
        atmosphere[1],  # N2
        atmosphere[2],  # O2
        atmosphere[6],  # Ar
        atmosphere[5],  # H
        atmosphere[7],  # N
        atmosphere[8],  # Anomalous O
        atmosphere[9],  # NO
        atmosphere[0]   # Total
    ])

    # Assign total kinetic temperature from model
    param_eq["Tinf"] = atmosphere[-1]
    param_eq["rho"] = rho

    # Calculate mean molecular mass [g/mol]
    total_density = np.sum(rho[:8])  # Exclude NO and total
    mmean = (
        rho[0] * constants.mHe + rho[1] * constants.mO + rho[2] * constants.mN2 +
        rho[3] * constants.mO2 + rho[4] * constants.mAr + rho[5] * constants.mH +
        rho[6] * constants.mN + rho[7] * constants.mAnO + rho[8] * constants.mNO
    ) / total_density

    param_eq["mmean"] = mmean
    param_eq["massConc"] = rho[:8] / np.sum(rho[:8])  # Mass fractions (first 8 species)

    # Specific gas constant [J/kg·K]
    param_eq["Rmean"] = (constants.R / mmean) * 1000

    # Orbital velocity assuming circular orbit [m/s]
    param_eq["vinf"] = np.sqrt(constants.mu_E / (constants.R_E + h))

    # Thermal speed [m/s]
    v_th = np.sqrt(2 * constants.kb * param_eq["Tinf"] / (mmean / constants.NA / 1000))
    param_eq["vth"] = v_th

    # Speed ratio: vinf / vth
    param_eq["s"] = param_eq["vinf"] / v_th

    return param_eq

