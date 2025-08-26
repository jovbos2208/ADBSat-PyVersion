import numpy as np
import pymsis
from .ADBSatConstants import ConstantsData

def get_atmospheric_data(param_eq, date, lon, lat, alt, f107, f107a, ap):
    """
    Calculates and enriches GSI model parameters with atmospheric data using pymsis (NRLMSIS2.1).

    This function calls the MSIS model for a specific time and location, then computes
    derived environmental parameters needed for aerodynamic calculations.

    Parameters:
        param_eq (dict): Base parameter dictionary to update.
        date (datetime): The date and time for the atmospheric calculation.
        lon (float): Geodetic longitude [degrees].
        lat (float): Geodetic latitude [degrees].
        alt (float): Altitude [km].
        f107 (float): Daily F10.7 solar flux index.
        f107a (float): 3-month average of F10.7 solar flux.
        ap (float or np.ndarray): Magnetic index (daily or an array of 7 values).

    Returns:
        dict: Updated `param_eq` dictionary with computed atmospheric values.
    """
    constants = ConstantsData()
    h = alt * 1000.0  # Convert altitude to meters for calculations

    # Create a list of 25 switches, all set to False (or 0) initially.
    # The default behavior of MSIS uses a value of 1 for all switches.
    # We will set switch 9 to 1 to enable NO calculation.
    options = [1] * 25
    options[9] = 1 # Switch for NO

    # Get atmospheric data from NRLMSIS2.1 model
    # The output is a 2D numpy array with shape (1, 11) for a single point calculation.
    # We access the first row for our data point, then columns by index.
    # Order: 0:Total, 1:He, 2:O, 3:N2, 4:O2, 5:Ar, 6:H, 7:N, 8:Anom. O, 9:NO, 10:T_exo
    msis_data = pymsis.calculate(date, lon, lat, alt, f107, f107a, ap, options=options)

    # Create a density vector consistent with the original 'environment.py' structure
    # Densities are in [m^-3], except for 'Total' which is mass density [kg/m^3]
    rho = np.array([
        msis_data[0, 1],  # He
        msis_data[0, 2],  # O
        msis_data[0, 3],  # N2
        msis_data[0, 4],  # O2
        msis_data[0, 5],  # Ar
        msis_data[0, 6],  # H
        msis_data[0, 7],  # N
        msis_data[0, 8],  # Anomalous O
        msis_data[0, 9],  # NO
        msis_data[0, 0]   # Total mass density
    ])

    # Assign kinetic temperature and densities to the parameter dictionary
    param_eq["Tinf"] = msis_data[0, 10]  # Exospheric temperature
    param_eq["rho"] = rho

    # Calculate mean molecular mass [g/mol]
    # Sum of (number_density_i * mass_i) / sum_of_number_densities
    number_densities = rho[:9]  # Species number densities
    total_number_density = np.sum(number_densities)
    
    molecular_masses = np.array([
        constants.mHe, constants.mO, constants.mN2, constants.mO2, 
        constants.mAr, constants.mH, constants.mN, constants.mAnO, constants.mNO
    ])
    
    mmean = np.sum(number_densities * molecular_masses) / total_number_density
    param_eq["mmean"] = mmean

    # Calculate mass concentrations for the first 8 species (used by DRIA model)
    mass_densities = number_densities[:8] * molecular_masses[:8]
    param_eq["massConc"] = mass_densities / np.sum(mass_densities)

    # Specific gas constant [J/(kg·K)]
    param_eq["Rmean"] = (constants.R / mmean) * 1000

    # Orbital velocity assuming a circular orbit [m/s]
    param_eq["vinf"] = np.sqrt(constants.mu_E / (constants.R_E + h))

    # Most probable thermal speed [m/s]
    mmean_kg = mmean / constants.NA / 1000  # Mean molecular mass in kg
    v_th = np.sqrt(2 * constants.kb * param_eq["Tinf"] / mmean_kg)
    param_eq["vth"] = v_th

    # Speed ratio (V_inf / v_th)
    param_eq["s"] = param_eq["vinf"] / v_th

    return param_eq
