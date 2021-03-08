import numpy as np
from scipy.stats import spearmanr, pearsonr

def clip(x):
    return np.clip(x, 1e-5, 1-1e-5)


def logit(x):
    x = clip(x)
    return np.log(x) - np.log(1 - x)

def na2minus1(a):
    inds = np.where(np.isnann(a))
    a[inds] = -1
    return a

# Define NA omit correlation metrics
# If none-NA value less than 2, return 1
def nanspearmanr(x, y):
    na = np.logical_or(np.isnan(x), np.isnan(y))
    x = x[~na]
    y = y[~na]
    if len(x) < 2:
        return (np.nan, np.nan)
    else:
        return spearmanr(x, y)
        
def nanpearsonr(x, y):
    x = x.flatten()
    y = y.flatten()
    na = np.logical_or(np.isnan(x), np.isnan(y))
    x = x[~na]
    y = y[~na]
    if len(x) < 2:
        return (np.nan, np.nan)
    else:
        return pearsonr(x, y)

def nanspearmanr_matrix(x, y):
    rho = np.zeros(x.shape[0], dtype=np.float)
    for i in range(x.shape[0]):
        na = np.logical_or(np.isnan(x[i]), np.isnan(y[i]))
        xi = x[i][~na]
        yi = y[i][~na]
        if len(xi) < 2:
            rho[i] = 1.
        else:
            rho[i] = spearmanr(xi, yi)[0]
    return rho
