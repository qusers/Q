import math
import numpy as np

from scipy import stats
from sklearn.utils import resample
from sklearn.metrics import mean_absolute_error
from sklearn.metrics import matthews_corrcoef

def pearsonr_ci(x, y, alpha=0.05):
    """ calculate Pearson correlation along with the confidence interval using scipy and numpy
    Parameters
    ----------
    x, y : iterable object such as a list or np.array
    Input for correlation calculation
    alpha : float
    Significance level. 0.05 by default
    Returns
    -------
    r : float
    Pearson's correlation coefficient
    pval : float
    The corresponding p value
    lo, hi : float
    The lower and upper bound of confidence intervals
    """
    x = np.array(x)
    y = np.array(y)
    r, p = stats.pearsonr(x, y)
    r_z = np.arctanh(r)
    se = 1 / math.sqrt(x.size - 3)
    z = stats.norm.ppf(1 - alpha / 2)
    lo_z, hi_z = r_z - z * se, r_z + z * se
    lo, hi = np.tanh((lo_z, hi_z))
    return r, p, lo, hi

def rmse(targets, predictions):
    """calculate root mean square error between array of targets and predictions"""
    return math.sqrt(((np.array(predictions) - np.array(targets)) ** 2).mean())

def analysis(X, Y, Z):
    """Returns statistics for correlation between arrays of targets and predictions"""
    n_size = int(len(X) * 0.90)
    bootstrap_dict = {'MCC': [], 'MAE': [], 'RMSE': [], 'k_tau': [], 'r2': [], 'lo': [], 'hi': [], 'p_value': [], 'SEM': []}
    data = np.column_stack((X, Y, Z))
    n = len(X)
    if n > 2:
        slope, intercept, _, p_value, std_err = stats.linregress(X, Y)
        regr = [(slope * x + intercept) for x in X]
        rho, _ = stats.spearmanr(X, Y)
        for i in range(1000):
            sample = resample(data, n_samples=n_size, random_state=i)
            X = [i[0] for i in sample]
            Y = [i[1] for i in sample]
            Z = [i[2] for i in sample]
            X2 = np.sign(X)
            Y2 = np.sign(Y)
            if all(i < 0 for i in X2) and all(j < 0 for j in Y2):
                MCC = 0.0
            else:
                bootstrap_dict['MCC'].append(matthews_corrcoef(X2, Y2))
            bootstrap_dict['MAE'].append(mean_absolute_error(X, Y))
            bootstrap_dict['RMSE'].append(rmse(X, Y))
            bootstrap_dict['k_tau'].append(stats.kendalltau(X, Y)[0])
            bootstrap_dict['p_value'].append(stats.kendalltau(X, Y)[1])
            r, _, lo, hi = pearsonr_ci(X, Y)
            bootstrap_dict['lo'].append(lo**2)
            bootstrap_dict['hi'].append(hi**2)
            bootstrap_dict['r2'].append(r**2)
            bootstrap_dict['SEM'].append(np.mean(Z))
        metrics = {'MCC': [], 'MAE': [], 'RMSE': [], 'sp_rho': rho, 'k_tau': [], 'r2': [], 'lo': [], 'hi': [], 'p_value': [], 'slope': slope, 'intercept': intercept, "n": n, "reg": regr, 'SEM': []}
        for k, v in bootstrap_dict.items():
            metrics[k].extend([round(np.mean(v), 2), round(np.std(v), 2)])
            #metrics[k].append(round(np.mean(v), 2))
            #metrics[k].append(round(np.std(v), 2))
        for k, v in metrics.items():
            try:
                metrics[k] = round(v, 2)
            except TypeError:
                continue
        return metrics