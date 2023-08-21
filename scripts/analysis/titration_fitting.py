import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import click
from scipy import optimize

sns.set_style("white")
sns.set_context("talk")


def curve_func_new(conc, K, n):
    return ((conc / K) ** n) / (1 + (conc / K) ** n)


def fit_bootstrap(p0, datax, datay, function, yerr_systematic=0.0):
    errfunc = lambda p, x, y: function(x, p[0], p[1]) - y
    # Fit first time
    pfit, perr = optimize.leastsq(errfunc, p0, args=(datax, datay), full_output=0)
    # Get the stdev of the residuals
    residuals = errfunc(pfit, datax, datay)
    sigma_res = np.std(residuals)
    sigma_err_total = np.sqrt(sigma_res**2 + yerr_systematic**2)
    # 100 random data sets are generated and fitted
    ps = []
    for i in range(100):
        randomDelta = np.random.normal(0.0, sigma_err_total, len(datay))
        randomdataY = datay + randomDelta
        randomfit, randomcov = optimize.leastsq(
            errfunc, p0, args=(datax, randomdataY), full_output=0
        )
        ps.append(randomfit)
    ps = np.array(ps)
    mean_pfit = np.mean(ps, 0)
    # You can choose the confidence interval that you want for your
    # parameter estimates:
    Nsigma = 1.0  # 1sigma gets approximately the same as methods above
    # 1sigma corresponds to 68.3% confidence interval
    # 2sigma corresponds to 95.44% confidence interval
    err_pfit = Nsigma * np.std(ps, 0)

    pfit_bootstrap = mean_pfit
    perr_bootstrap = err_pfit
    return pfit_bootstrap, perr_bootstrap


def fit_titration(x, y):
    p0 = [1, 1]
    pfit, perr = fit_bootstrap(p0, x, y, curve_func_new, yerr_systematic=0.0)
    return pfit, perr


def normalize_data_full(data):
    if np.min(data) == np.max(data):
        return data
    return (data - np.min(data)) / (np.max(data) - np.min(data))


def normalize_data(data):
    if np.min(data) == np.max(data):
        return data
    return (data) / (np.max(data))


def compute_mg_1_2(mg_conc, data):
    pstart = [1, 1]
    norm_data = -normalize_data(np.array(data)) + 1
    pfit, perr = fit_bootstrap(pstart, mg_conc, norm_data, curve_func_new)
    return pfit, perr


@click.command()
@click.argument("json_file", type=click.Path(exists=True))
def main(json_file):
    """
    main function for script
    """
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
