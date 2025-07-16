import numpy as np
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
    Routes computation to the correct GSI coefficient model.

    Parameters:
        param_eq (dict): Model parameters
        delta (np.ndarray): Panel angles to freestream [rad]
        matID (np.ndarray): Material identifiers (unused here)

    Returns:
        Tuple[np.ndarray]: cp, ctau, cd, cl per panel
    """
    model = param_eq.get('gsi_model', None)
    if model is None:
        raise ValueError("Missing 'gsi_model' in parameters.")

    # Dispatch to model
    model = model.lower()
    if model == "cll":
        return coeff_cll(param_eq, delta)
    elif model == "cook":
        return coeff_cook(param_eq, delta)
    elif model == "dria":
        return coeff_DRIA(param_eq, delta)
    elif model == "maxwell":
        return coeff_maxwell(param_eq, delta)
    elif model == "newton":
        return coeff_newton(param_eq, delta)
    elif model == "schaaf":
        return coeff_schaaf(param_eq, delta)
    elif model == "sentman":
        return coeff_sentman(param_eq, delta)
    elif model == "storchhyp":
        return coeff_storchHyp(param_eq, delta)
    else:
        raise ValueError(f"Unknown GSI model: {model}")

