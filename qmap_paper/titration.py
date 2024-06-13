"""
Handles fitting of 16 point mg2+ titrations for estimating the Mg2+/2 using 
the hill coefficient equation.
"""

import numpy as np
from scipy import optimize
from typing import List, Tuple


def normalized_hill_equation(conc: float, K: float, n: float, A: float) -> float:
    """
    Calculate the normalized Hill equation value.

    Args:
        conc (float): Concentration of the titration agent (mg2+ or a ligand).
        K (float): Dissociation constant.
        n (float): Hill coefficient.
        A (float): Maximum value.

    Returns:
        float: The calculated normalized Hill equation value.
    """
    return A * ((conc / K) ** n) / (1 + (conc / K) ** n)


def fit_bootstrap(
    p0: List[float],
    x: np.ndarray,
    y: np.ndarray,
    function: callable,
    n_runs: int = 100,
    n_sigma: float = 1.0,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Uses bootstrap method to estimate the 1 sigma confidence interval of the parameters
    for a fit to data (x,y) with function function(x, params).

    1 sigma corresponds to 68.3% confidence interval
    2 sigma corresponds to 95.44% confidence interval

    Args:
        p0: Initial guess for parameters.
        x: Independent values.
        y: Dependent values, what you are trying to fit to.
        function: Function to fit to, should be a python function.
        n_runs: Number of bootstrap runs. Default is 100.
        n_sigma: Number of sigma to use for confidence interval. Default is 1.0.

    Returns:
        Tuple containing the mean parameter values and the parameter errors.

    Raises:
        None

    Examples:
        p0 = [1.0, 2.0, 3.0]
        x = np.array([1, 2, 3, 4, 5])
        y = np.array([2, 4, 6, 8, 10])
        function = lambda x, a, b, c: a*x**2 + b*x + c
        mean_pfit, err_pfit = fit_bootstrap(p0, x, y, function)
    """
    errfunc = lambda p, x, y: function(x, p[0], p[1], p[2]) - y
    # Fit first time
    pfit, perr = optimize.leastsq(errfunc, p0, args=(x, y), full_output=0)
    # Get the stdev of the residuals
    residuals = errfunc(pfit, x, y)
    sigma_res = np.std(residuals)
    sigma_err_total = np.sqrt(sigma_res**2)
    # 100 random data sets are generated and fitted
    ps = []
    for i in range(n_runs):
        random_delta = np.random.normal(0.0, sigma_err_total, len(y))
        random_y = y + random_delta
        random_fit, _ = optimize.leastsq(errfunc, p0, args=(x, random_y))
        # a hack to stop insane values
        if random_fit[2] > 10:
            continue
        ps.append(random_fit)
    ps = np.array(ps)
    mean_pfit = np.mean(ps, 0)
    err_pfit = n_sigma * np.std(ps, 0)
    return mean_pfit, err_pfit


def normalize_data_full(data: np.ndarray) -> np.ndarray:
    """
    Normalize the given data using the full range of values.

    Args:
        data: A numpy array containing the data to be normalized.

    Returns:
        A numpy array with the normalized data.

    Raises:
        None

    Examples:
        >>> data = np.array([1, 2, 3, 4, 5])
        >>> normalize_data_full(data)
        array([0.  , 0.25, 0.5 , 0.75, 1.  ])
    """
    if np.min(data) == np.max(data):
        return data
    return (data - np.min(data)) / (np.max(data) - np.min(data))


def normalize_data(data: np.ndarray) -> np.ndarray:
    """
    Normalize the input data array.

    Args:
        data: A numpy array containing the data to be normalized.

    Returns:
        A numpy array with the normalized data.

    Raises:
        None.

    Examples:
        >>> data = np.array([1, 2, 3, 4, 5])
        >>> normalize_data(data)
        array([0. , 0.25, 0.5 , 0.75, 1. ])
    """
    if np.min(data) == np.max(data):
        return data
    return data / np.max(data)


def compute_mg_1_2(
    mg_conc: List[float], data: List[float]
) -> Tuple[List[float], List[float]]:
    """
    Computes the fitting parameters for the normalized Hill equation.

    Args:
        mg_conc (List[float]): A list of magnesium concentrations.
        data (List[float]): A list of experimental data.

    Returns:
        Tuple[List[float], List[float]]: A tuple containing the fitting parameters (pfit) and their errors (perr).
    """
    pstart = [1, 1, 0.5]
    norm_data = -normalize_data(np.array(data)) + 1
    pfit, perr = fit_bootstrap(pstart, mg_conc, norm_data, normalized_hill_equation)
    return pfit, perr
