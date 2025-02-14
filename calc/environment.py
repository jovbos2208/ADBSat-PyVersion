import numpy as np
from .ADBSatConstants import ConstantsData,EnvironmentData
import pandas as pd



def environment(param_eq,database,idx,h):
    """
    Computes a number of atmospheric parameters assuming a circular orbit without atmospheric co-rotation or winds.

    Parameters:
        h (float): Altitude [m]
        lat (float): Geodetic latitude [deg]
        lon (float): Longitude [deg]
        day_of_year (int): Day of the year (1-365)
        ut_seconds (float): Seconds of the day
        f107_avg (float): 81-day average of F10.7 flux (centered on day_of_year)
        f107_daily (float): Daily F10.7 flux for the previous day
        magnetic_index (list): Magnetic index (AP) information [1x7 array]
        ano_flag (bool): Flag for anomalous oxygen (True/False)

    Returns:
        dict: Dictionary containing atmospheric parameters.
    """
    constants = ConstantsData()

    # Convert altitude to km for MSISE-00
    alt_km = h / 1000.0


    

    #database_idx = int(np.random.uniform(0,database.shape[0]))
    
    #print(database_idx)

    # Run MSISE-00 atmospheric model
    atmosphere = database.iloc[idx].to_numpy()

    # Extract outputs from the model
    # Extract exospheric and kinetic temperatures
    rho = np.array([
        atmosphere[4], atmosphere[3], atmosphere[1],
        atmosphere[2], atmosphere[6], atmosphere[5],
        atmosphere[7], atmosphere[8], atmosphere[9], atmosphere[0]  # Adjust for other species if needed
    ])
    #print(rho)
    param_eq["Texo"] = atmosphere[-2]
    param_eq["Tinf"] = atmosphere[-1]
    param_eq["rho"] = rho
    

    # Calculate mean molecular mass [g/mol]
    
    total_density = np.sum(rho[:8]) 
    mmean = (rho[0] * constants.mHe + rho[1] * constants.mO +
             rho[2] * constants.mN2 + rho[3] * constants.mO2 +
             rho[4] * constants.mAr + rho[5] * constants.mH +
             rho[6] * constants.mN + rho[7] * constants.mAnO + rho[8] * constants.mNO) / total_density

    param_eq["mmean"] = mmean
    

    # Calculate mass concentrations
    param_eq["massConc"] = np.zeros(10)
    param_eq["massConc"][0] = rho[0] / rho[9] * (constants.mHe / constants.NA / 1000)
    param_eq["massConc"][1] = rho[1] / rho[9] * (constants.mO / constants.NA / 1000)
    param_eq["massConc"][2] = rho[2] / rho[9] * (constants.mN2 / constants.NA / 1000)
    param_eq["massConc"][3] = rho[3] / rho[9] * (constants.mO2 / constants.NA / 1000)
    param_eq["massConc"][4] = rho[4] / rho[9] * (constants.mAr / constants.NA / 1000)
    param_eq["massConc"][5] = rho[5] / rho[9] * (constants.mH / constants.NA / 1000)
    param_eq["massConc"][6] = rho[6] / rho[9] * (constants.mN / constants.NA / 1000)
    param_eq["massConc"][7] = rho[7] / rho[9] * (constants.mN / constants.NA / 1000)
    param_eq["massConc"][8] = rho[8] / rho[9] * (constants.mNO / constants.NA / 1000)

    # Calculate specific gas constant [J/kg K]
    param_eq["Rmean"] = (constants.R / mmean) * 1000

    # Orbital velocity [m/s]
    param_eq["vinf"] = np.sqrt(constants.mu_E / (constants.R_E + h))

    # Thermal velocity
    param_eq["vth"] = np.sqrt(2 * constants.kb * param_eq["Tinf"] / (mmean / constants.NA / 1000))

    # Speed ratio
    param_eq["s"] = param_eq["vinf"] / param_eq["vth"]

    return param_eq
