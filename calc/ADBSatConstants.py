import numpy as np

class ConstantsData:
    """
    Container for universal physical and gas-specific constants.
    Units are consistent with SI unless otherwise noted.
    """
    def __init__(self):
        # Earth and physical constants
        self.mu_E = 3.986004418e14  # Gravitational parameter (Earth) [m^3/s^2]
        self.R_E = 6.37813649e6     # Equatorial radius of Earth [m]
        self.R = 8.31446261815324   # Universal gas constant [J/(mol·K)]
        self.kb = 1.3806503e-23     # Boltzmann constant [J/K]
        self.NA = 6.02214076e23     # Avogadro's number [1/mol]
        self.pi = np.pi             # Mathematical constant π

        # Molecular masses [g/mol]
        self.mHe = 4.002602
        self.mO = 15.9994
        self.mN2 = 28.0134
        self.mO2 = 31.9988
        self.mAr = 39.948
        self.mH = 1.0079
        self.mN = 14.0067
        self.mAnO = 15.99       # Anomalous oxygen
        self.mNO = 30.0067      # Nitric oxide


class EnvironmentData:
    """
    Structure for storing environment-related input parameters.
    Most are expected to be scalar floats or NumPy arrays.
    """
    def __init__(self):
        self.alpha = np.array([])       # Energy accommodation coefficient
        self.alphaN = np.array([])      # Normal accommodation coefficient
        self.sigmaN = np.array([])      # Normal momentum accommodation coefficient
        self.sigmaT = np.array([])      # Tangential momentum accommodation coefficient
        self.Tw = np.array([])          # Wall temperature [K]
        self.Tinf = np.array([])        # Kinetic temperature [K]
        self.vinf = np.array([])        # Incident velocity [m/s]
        self.Rmean = np.array([])       # Specific gas constant [J/(kg·K)]
        self.s = np.array([])           # Speed ratio (V_inf / v_th)
        self.Vw = np.array([])          # Mean thermal velocity of reflected particles [m/s]
        self.gamma = np.array([])       # cos(δ), used in DRIA model
        self.ell = np.array([])         # sin(δ), used in DRIA model
        self.massConc = np.array([[]])  # Species mass concentrations (8×N panels)
        self.sol_cR = np.array([[]])    # Specular reflectivity for SRP model
        self.sol_cD = np.array([[]])    # Diffuse reflectivity for SRP model
        self.gsi_model = ""             # GSI model identifier (string)

