import numpy as np
import time


def periods_statistic(evnts, intervals, bins, pmin, pmax, n_steps=1000, exposure_on=True):

    def _chi_2(fold, evnts, intervals, exposure, periods, bins):
        print('start folding')
        folded = _fold(evnts, periods, bins)  # (n_periods, bins)
        print('folded')
        print('start exposuring')
        expo = exposure(intervals, periods, bins)
        print('exposured')
        folded1 = folded/expo
        folded2 = folded*(2-expo)

        expected = np.sum(folded, axis=1)/bins  # (n_periods, 1)
        chi_square = np.sum((folded.T - expected)**2/expected, axis=0)

        expected1 = np.sum(folded1, axis=1)/bins  # (n_periods, 1)
        chi_square1 = np.sum((folded1.T - expected1)**2/expected1, axis=0)

        expected2 = np.sum(folded2, axis=1)/bins  # (n_periods, 1)
        chi_square2 = np.sum((folded2.T - expected2)**2/expected2, axis=0)

        return chi_square.T, chi_square1.T, chi_square2.T

    def _fold(events, periods, bins):
        n = events.shape[0]
        m = periods.shape[0]
        e = np.repeat(events, m).reshape(n, m)  # (events, periods)
        x = np.mod(e, periods)/periods
        folded = np.array([np.count_nonzero(x//(1/bins) == i, axis=0)
                           for i in range(bins)])
        return folded.T

    def _exposure(intervals, periods, bins):
        res = np.array([])

        for period in periods:
            dt = period/bins
            exp_times = np.zeros(bins)

            i = (np.ceil(intervals.transpose()[0]/dt) % bins).astype(int)
            indx = 0
            for i1 in i:
                exp_times[i1] += dt * \
                    np.ceil(intervals[indx][0]/dt) - intervals[indx][0]
                indx += 1

            j = (np.floor(intervals.transpose()[1]/dt) % bins).astype(int)
            indx = 0
            for j1 in j:
                exp_times[j1] += intervals[indx][1] - \
                    dt*np.floor(intervals[indx][1]/dt)
                indx += 1

            num_periods = np.count_nonzero(
                ((dt*np.floor(intervals.transpose()[1]/dt) - dt*np.ceil(intervals.transpose()[0]/dt))//period) >= 0)
            exp_times += num_periods*dt

            for k in range(len(i)):
                if i[k] < j[k]:
                    exp_times[i[k]: j[k]] += dt

                elif i[k] > j[k]:
                    exp_times[i[k]:] += dt
                    exp_times[:j[k]] += dt

            exp_times /= np.max(exp_times)
            res = np.concatenate((res, exp_times))
        return res.reshape((periods.shape[0], bins))

    start_time = time.time()
    periods = np.linspace(pmin, pmax, n_steps)
    stat, stat1, stat3 = _chi_2(
        _fold, evnts, intervals, _exposure, periods, bins)  # (n_periods, 1)
    print("Completed for ", (time.time()-start_time)/60, " minutes")
    return periods, stat, stat1, stat3
