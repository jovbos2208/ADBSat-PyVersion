import numpy as np

class ConstantsData:
    """
    Class to define universal constants.
    """
    def __init__(self):
        self.mu_E = 3.986004418e14  # GM Earth [m^3/s^2]
        self.R_E = 6.37813649e6     # Equatorial Earth radius [m]
        self.R = 8.31446261815324   # Universal Gas constant [J K^-1 mol^-1]
        self.kb = 1.3806503e-23     # Boltzmann constant [m^2 kg^-2 K^-1]
        self.NA = 6.02214076e23     # Avogadro constant [n mol^-1]
        self.pi = np.pi             # Pi
        
        # Molecular masses [g/mol]
        self.mHe = 4.002602
        self.mO = 15.9994
        self.mN2 = 28.0134
        self.mO2 = 31.9988
        self.mAr = 39.948
        self.mH = 1.0079
        self.mN = 14.0067
        self.mAnO = 15.99
        self.mNO = 30.0067


class EnvironmentData:
    """
    Class to define environmental data structure.
    """
    def __init__(self):
        # Initialize all parameters as empty arrays
        self.alpha = np.array([])       # Energy accommodation coefficient
        self.alphaN = np.array([])      # Normal accommodation coefficient
        self.sigmaN = np.array([])      # Normal accommodation coefficient [-]
        self.sigmaT = np.array([])      # Tangential accommodation coefficient [-]
        self.Tw = np.array([])          # Wall temperature [K]
        self.Tinf = np.array([])        # Kinetic temperature [K]
        self.vinf = np.array([])        # Incident velocity [m/s]
        self.Rmean = np.array([])       # Specific gas constant [J/kg K]
        self.s = np.array([])           # Thermal speed ratio [-]
        self.Vw = np.array([])          # Average normal velocity of diffusely reflected molecules [m/s]
        self.gamma = np.array([])       # Auxiliary variable for DRIA model
        self.ell = np.array([])         # Auxiliary variable for DRIA model
        self.massConc = np.array([[]])  # Mass concentration of species
        self.sol_cR = np.array([[]])
        self.sol_cD = np.array([[]])
        self.gsi_model = ""             # GSI model name
