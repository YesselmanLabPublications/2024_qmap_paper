"""
Handles fitting of 16 point mg2+ titrations for estimating the Mg2+/2 using 
the hill coefficient equation.
"""

import numpy as np
from scipy import optimize


def normalized_hill_equation(conc, K, n, A):
    """
    Assumes a range from 0 to 1
    :param conc: concentration of titration agent either mg2+ or a ligand
    :param K: dissociation constant
    :param n: hill coefficient
    :param A: maximum value
    """
    return A * ((conc / K) ** n) / (1 + (conc / K) ** n)


def fit_bootstrap(p0, x, y, function, n_runs=100, n_sigma=1.0):
    """
    Uses bootstrap method to estimate the 1 sigma confidence interval of the parameters
    for a fit to data (x,y) with function function(x, params).

    1 sigma corresponds to 68.3% confidence interval
    2 sigma corresponds to 95.44% confidence interval

    :param p0: initial guess for parameters
    :param x: x - independent values
    :param y: y - dependent values, what are you trying to fit to
    :param function: function to fit to, should be a python function
    :param n_runs: number of bootstrap runs - default 100
    :param n_sigma: number of sigma to use for confidence interval - default 1.0
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
        random_fit, _ = optimize.leastsq(errfunc, p0, args=(x, random_y), full_output=0)
        ps.append(random_fit)
    ps = np.array(ps)
    mean_pfit = np.mean(ps, 0)
    err_pfit = n_sigma * np.std(ps, 0)
    return mean_pfit, err_pfit


def normalize_data_full(data):
    if np.min(data) == np.max(data):
        return data
    return (data - np.min(data)) / (np.max(data) - np.min(data))


def normalize_data(data):
    if np.min(data) == np.max(data):
        return data
    return (data) / (np.max(data))


def compute_mg_1_2(mg_conc, data):
    pstart = [1, 1, 1]
    norm_data = -normalize_data(np.array(data)) + 1
    pfit, perr = fit_bootstrap(pstart, mg_conc, norm_data, normalized_hill_equation)
    return pfit, perr


"""
@click.command()
@click.argument("json_file", type=click.Path(exists=True))
def main(json_file):
    df = pd.read_json(json_file)
    df = df[df["mg_conc"] != 5.0]
    data = []
    for name, g in df.groupby("name"):
        r = compute_mg_1_2(list(g["mg_conc"]), list(g["gaaa_avg"]))
        data.append([name, r[0][0], r[0][1]])
    df_r = pd.DataFrame(data, columns="name,k,n".split(","))
    df_r.to_csv("data/mg_1_2.csv", index=False)


# pylint: disable=no-value-for-parameter
if __name__ == "__main__":
    main()
"""
