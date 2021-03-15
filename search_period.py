import numpy as np
import time


def periods_statistic(events, intervals, nbins, pmin, pmax, n_steps=1000):

    def _chi_2(fold, evnts, intervals, exposure, periods, nbins):
        folded = _fold(evnts, periods, nbins)  # (n_periods, nbins)

        if intervals.shape[0] == 1:
            expected = np.sum(folded, axis=1)/nbins  # (n_periods, 1)
            chi_square = np.sum((folded.T - expected)**2/expected, axis=0)
            return chi_square, chi_square

        expected = np.sum(folded, axis=1)/nbins  # (n_periods, 1)
        chi_square = np.sum((folded.T - expected)**2/expected, axis=0)

        expo = exposure(intervals, periods, nbins)  # (n_periods, nbins)
        # the average burst rate per unit exposure time
        p = len(evnts)/np.sum(intervals.T[1]-intervals.T[0])
        expected_expo = p*expo
        chi_square_expo = np.zeros(len(periods))
        for i, period in enumerate(periods):
            for j in range(nbins):
                if expo[i][j] > 0:
                    chi_square_expo[i] += (folded[i][j] -
                                           expected_expo[i][j])**2/expected_expo[i][j]
        # expo_time = exposure(intervals, periods, nbins)  # (n_periods, nbins)
        # expo = expo = (expo_time.T/expo_time.max(axis=1)).T
        # expected_expo = len(evnts)/np.sum(intervals.T[1]-intervals.T[0])*expo
        # chi_square_expo = np.zeros(len(periods))

        # expected1 = np.sum(folded1, axis=1)/nbins  # (n_periods, 1)
        # chi_square1 = np.sum((folded1.T - expected)**2/expected1, axis=0)

        # expected2 = len(evnts)/np.sum(intervals.T[1]-intervals.T[0])*expo
        # chi_square2 = np.sum((folded2 - expected2)**2/expected2, axis=1)
        return chi_square, chi_square_expo

    def _fold(events, periods, nbins):
        folded = np.zeros((len(periods), nbins))
        for i, period in enumerate(periods):
            folded[i] = np.histogram(
                events % period/period, nbins, range=(0, 1))[0]
        return folded

    def _exposure(intervals, periods, nbins):
        res = np.zeros((periods.shape[0], nbins))

        for indx, period in enumerate(periods):
            dt = period/nbins
            exp_times = np.zeros(nbins)

            for interval in intervals:
                i = int(np.floor(interval[0]/dt) % nbins)
                t = dt*np.ceil(interval[0]/dt)
                if t > interval[1]:
                    exp_times[i] += interval[1]-interval[0]
                    continue
                elif t+dt >= interval[1]:
                    exp_times[i] += t-interval[0]
                    exp_times[(i+1) % nbins] += interval[1] - t
                    continue
                else:
                    exp_times[i] += t-interval[0]
                while t <= interval[1]:
                    i = (i+1) % nbins
                    if t+dt < interval[1]:
                        exp_times[i] += dt
                    else:
                        exp_times[i] += interval[1] - t
                    t += dt
            res[indx] += exp_times
        return res.reshape((periods.shape[0], nbins))

    start_time = time.time()
    periods = np.linspace(pmin, pmax, n_steps)
    stat, stat_expo = _chi_2(
        _fold, events, intervals, _exposure, periods, nbins)  # (n_periods, 1)
    print("Completed for", (time.time()-start_time)/60, "minutes")
    return periods, stat, stat_expo
    # stat - без учета экспозиции, stat_expo - с учетом
