import numpy as np


def fold(events, period, bins):
    x = np.sort(events % period/period)
    folded = np.zeros(bins, dtype=int)
    for i in x:
        folded[int(i//(1/bins))] += 1
    return folded


def chi_2(evnts, fold, period, bins):
    #     expected = (len(evnts)/(evnts[-1]-evnts[0]))*period/bins
    folded = fold(evnts, period, bins)
    expected = np.sum(folded)/bins
    chi_square = 0
    for f in folded:
        chi_square += (f-expected)**2/expected
    return chi_square


def periods_statistic(evnts, chi_2, fold, bins, pmin, pmax):
    T = evnts[-1] - evnts[0]
    fmin = 1/pmax
    fmax = 1/pmin
    freq = fmin
    stat = []
    periods = []
    while freq < fmax:
        stat.append(chi_2(evnts, fold, 1/freq, bins))
        periods.append(1/freq)
        freq += 0.1/T
    return periods, stat
