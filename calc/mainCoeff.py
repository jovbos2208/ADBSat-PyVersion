import numpy as np
from .ADBSatConstants import ConstantsData
from .coeff_CLL import coeff_cll
from .coeff_cook import coeff_cook
from .coeff_DRIA import coeff_DRIA
from .coeff_maxwell import coeff_maxwell
from .coeff_newton import coeff_newton
from .coeff_schaaf import coeff_schaaf
from .coeff_sentman import coeff_sentman
from .coeff_storchHyp import coeff_storchHyp

def mainCoeff(param_eq, delta, matID):
    """
    Passes the input parameters to the correct GSI model.

    Parameters:
        param_eq (dict): Dictionary containing the parameters for the GSI model.
        delta (numpy.ndarray): Array of angles between the surface normal and the flow (radians).
        matID (numpy.ndarray): Array of material identifiers for each panel.

    Returns:
        tuple: cp, ctau, cd, cl (numpy arrays of coefficients for each panel).
    """
    # Determine the model to use
    gsi_model = param_eq.get('gsi_model', None)
    
    if gsi_model is None:
        raise ValueError("GSI model not specified in param_eq.")

    # Match the GSI model and compute coefficients
    if gsi_model == "CLL":
        cp, ctau, cd, cl = coeff_cll(param_eq, delta)
    elif gsi_model == "Cook":
        cp, ctau, cd, cl = coeff_cook(param_eq, delta)
    elif gsi_model == "DRIA":
        cp, ctau, cd, cl = coeff_DRIA(param_eq, delta)
    elif gsi_model == "Maxwell":
        cp, ctau, cd, cl = coeff_maxwell(param_eq, delta)
    elif gsi_model == "Newton":
        cp, ctau, cd, cl = coeff_newton(param_eq, delta)
    elif gsi_model == "Schaaf":
        cp, ctau, cd, cl = coeff_schaaf(param_eq, delta)
    elif gsi_model == "Sentman":
        cp, ctau, cd, cl = coeff_sentman(param_eq, delta)
    elif gsi_model == "StorchHyp":
        cp, ctau, cd, cl = coeff_storchHyp(param_eq, delta)
    else:
        raise ValueError(f"Unknown GSI model: {gsi_model}")

    return cp, ctau, cd, cl
