import numpy as np


def periods_statistic1(evnts, intervals, chi_2, fold, bins, pmin, pmax, n_steps=1000):
    freq = np.linspace(1/pmax, 1/pmin, n_steps)
    periods = 1/freq[::-1]
    stat = chi_2(evnts, fold, periods, bins)  # (n_periods, 1)
    return periods, stat


def chi_2(evnts, fold, periods, bins):
    folded = fold(evnts, periods, bins)  # (n_periods, bins)
    expected = np.sum(folded, axis=1)/bins  # (n_periods, 1)
    chi_square = np.sum((folded.T - expected)**2/expected, axis=0)
    return chi_square.T


def fold(events, periods, bins):
    n = events.shape[0]
    m = periods.shape[0]
    p = np.repeat(periods, n).reshape(m, n)
    e = np.repeat(events, m).reshape(n, m).T
    x = np.sort(e % p/p)
    folded = np.array([np.count_nonzero(x//(1/bins) == i, axis=1)
                       for i in range(bins)])
    return folded.T
