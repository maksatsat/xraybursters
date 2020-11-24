import numpy as np


def periods_statistic(evnts, intervals, bins, pmin, pmax, n_steps=1000, exposure_on=True):
    def chi_2(fold, evnts, exposure, periods, bins):
        folded = _fold(evnts, periods, bins)  # (n_periods, bins)
        if exposure_on:
            expo = exposure(intervals, periods, bins)
            folded = folded*(2-expo)

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

    def exposure(intervals, periods, bins):
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

            num_periods = (dt*np.floor(intervals.transpose()
                                       [1]/dt) - dt*np.ceil(intervals.transpose()[0]/dt))//period
            exp_times += np.sum(num_periods)*dt

            for k in range(len(i)):
                if i[k] < j[k]:
                    exp_times[i[k]: j[k]] += dt

                elif i[k] > j[k]:
                    exp_times[i[k]:] += dt
                    exp_times[:j[k]] += dt

            exp_times /= np.max(exp_times)
            res = np.concatenate((res, exp_times))
        return res.reshape((periods.shape[0], bins))
#     freq = np.linspace(1/pmax, 1/pmin, n_steps)
#     periods = 1/freq[::-1]
    periods = np.linspace(pmin, pmax, n_steps)
#     expo = _exposure(intervals, periods, bins)
    stat = chi_2(_fold, evnts, exposure, periods, bins)  # (n_periods, 1)
    return periods, stat
