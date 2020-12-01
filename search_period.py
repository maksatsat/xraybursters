import numpy as np
import time


def periods_statistic(evnts, intervals, bins, pmin, pmax, n_steps=1000, exposure_on=True):

    def _chi_2(fold, evnts, intervals, exposure, periods, bins):
        # print('start folding')
        folded = _fold(evnts, periods, bins)  # (n_periods, bins)
        # print('folded')
        # print('start exposuring')
        expo = exposure(intervals, periods, bins)  # (n_periods, bins)
        # print('exposured')
        folded1 = folded/((expo.T/np.max(expo, axis=1)).T)
        folded2 = folded*1

        expected = np.sum(folded, axis=1)/bins  # (n_periods, 1)
        chi_square = np.sum((folded.T - expected)**2/expected, axis=0)

        expected1 = np.sum(folded1, axis=1)/bins  # (n_periods, 1)
        chi_square1 = np.sum((folded1.T - expected)**2/expected1, axis=0)

        expected2 = len(evnts)/np.sum(intervals.T[1]-intervals.T[0])*expo
        chi_square2 = np.sum((folded2 - expected2)**2/expected2, axis=1)

        return chi_square.T, chi_square1.T, chi_square2.T

    def _fold(events, periods, bins):
        folded = np.zeros((len(periods), bins))
        for i in range(len(periods)):
            folded[i] = np.histogram(
                events % periods[i]/periods[i], bins, range=(0, 1))[0]
        return folded

    def _exposure_old(intervals, periods, bins):
        res = np.zeros((len(periods), bins))
        for i in range(len(periods)):
            for interval in intervals:
                res[i] += np.histogram(np.arange(
                    interval[0], interval[1], 1) % periods[i], bins, range=(0, periods[i]))[0]
        return res

    def _exposure(intervals, periods, bins):
        res = np.zeros((periods.shape[0], bins))

        indx = 0
        for period in periods:
            dt = period/bins
            exp_times = np.zeros(bins)

            for interval in intervals:
                i = int(np.floor(interval[0]/dt) % bins)
                t = dt*np.ceil(interval[0]/dt)
                if t > interval[1]:
                    exp_times[i] += interval[1]-interval[0]
                    continue
                elif t+dt >= interval[1]:
                    exp_times[i] += t-interval[0]
                    exp_times[(i+1) % bins] += interval[1] - t
                    continue
                else:
                    exp_times[i] += t-interval[0]
                while t <= interval[1]:
                    i = (i+1) % bins
                    if t+dt < interval[1]:
                        exp_times[i] += dt
                    else:
                        exp_times[i] += interval[1] - t
                    t += dt
            res[indx] += exp_times
            indx += 1
        return res.reshape((periods.shape[0], bins))

    start_time = time.time()
    periods = np.linspace(pmin, pmax, n_steps)
    stat, stat1, stat2 = _chi_2(
        _fold, evnts, intervals, _exposure, periods, bins)  # (n_periods, 1)
    print("Completed for", (time.time()-start_time)/60, "minutes")
    return periods, stat, stat1, stat2
