import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import datetime
from pathlib import Path
from astropy.time import Time
import pandas as pd
import interval
from event_generator import Events
from search_period import periods_statistic


def efsearch(name, nbins=10, n_steps=1000, num_simulations=100, saving=True, pmin=None, pmax=None):
    # saving_directory = 'results/' + name + '/'
    saving_directory = '/'
    bursts_url = 'https://burst.sci.monash.edu/aqoutput?dtype=bursts&fields=name%2Ctime%2Cinstr&output=topcat&timef=mjd&qfield=name&query_op=%3D&query_val=' + \
        name.replace(' ', '+')
    osbervations_url = 'https://burst.sci.monash.edu/aqoutput?dtype=observations&fields=name%2Ctstart%2Cinstr%2Ctstop&output=topcat&timef=mjd&qfield=name&query_op=%3D&query_val=' + \
        name.replace(' ', '+')

    print('Reading data: ', end='')
    bursts_table = pd.read_csv(
        bursts_url, comment='#', sep='\t', header=0, index_col=False)
    bursts_table['time_s'] = [datetime.datetime.strptime(
        x.iso, '%Y-%m-%d %H:%M:%S.%f') for x in Time(bursts_table['time'], format='mjd')]
    osbervations_table = pd.read_csv(
        osbervations_url, comment='#', sep='\t', header=0, index_col=False)
    osbervations_table['tstart_s'] = [datetime.datetime.strptime(
        x.iso, '%Y-%m-%d %H:%M:%S.%f') for x in Time(osbervations_table['tstart'], format='mjd')]
    osbervations_table['tstop_s'] = [datetime.datetime.strptime(
        x.iso, '%Y-%m-%d %H:%M:%S.%f') for x in Time(osbervations_table['tstop'], format='mjd')]
    instruments = bursts_table['instr'].unique()
    print('done')

    if saving:
        print('Saving: ', end='')
        burts_save = {
            'Bursts': bursts_table['time_s']
        }
        df_bursts = pd.DataFrame(burts_save)
        df_bursts.to_csv(saving_directory + name + ':bursts.csv', index=False)
        observations_save = {
            'Start': osbervations_table['tstart_s'],
            'Stop': osbervations_table['tstop_s']
        }
        df_observations = pd.DataFrame(burts_save)
        df_observations.to_csv(saving_directory + name +
                               ':observations.csv', index=False)
        print('done')

    print('Merging observation intervals: ', end='')
    # объединенные интервалы наблюдений
    merged_intervals = interval.interval()
    for i in range(len(instruments)):
        ob = osbervations_table[osbervations_table['instr'] == instruments[i]]
        ob.index = np.arange(len(ob))
        if (len(ob) > 0):
            for j in range(len(ob)):
                merged_intervals = merged_intervals | interval.interval(
                    [ob['tstart'][j], ob['tstop'][j]])
    print('done')
    print('Merging burst times: ', end='')
    # объединенные времена вспышек

    def notfoundsim(b, merged_bursts, delta_t=30):
        # if len(merged_bursts)==0:
            # return True
        for x in merged_bursts:
            delta = (Time(x, format='mjd')-Time(b, format='mjd')
                     ).to_datetime().total_seconds()
            if (np.abs(delta) <= delta_t):
                return False
        return True

    merged_bursts = []
    for i in range(len(instruments)):
        ob = osbervations_table[osbervations_table['instr'] == instruments[i]]
        if len(ob) > 0:
            b = bursts_table[bursts_table['instr'] == instruments[i]]['time']
            for j in b:
                if notfoundsim(j, merged_bursts):
                    merged_bursts.append(j)
    merged_bursts.sort()
    # Converting from MJD to Y-m-d h:m:s.f string
    time_intervals = []
    for interv in merged_intervals:
        tm1 = Time(interv[0], format='mjd')
        d1 = datetime.datetime.strptime(
            tm1.iso, '%Y-%m-%d %H:%M:%S.%f').strftime("%Y-%m-%d %H:%M:%S.%f")
        tm2 = Time(interv[1], format='mjd')
        d2 = datetime.datetime.strptime(
            tm2.iso, '%Y-%m-%d %H:%M:%S.%f').strftime("%Y-%m-%d %H:%M:%S.%f")
        time_intervals.append([d1, d2])
    # Converting from MJD to datetime.datetime
    events = []
    for burst in merged_bursts:
        tm = Time(burst, format='mjd')
        d = datetime.datetime.strptime(tm.iso, '%Y-%m-%d %H:%M:%S.%f')
        events.append(d)

    # start datetime.datetime - left interval
    # intrvls - intervals relatively to start in seconds
    # evnts - events relatively to start in seconds
    start, intrvls = Events.dates_to_seconds(time_intervals)
    evnts = []
    for event in events:
        evnts.append((event-start).total_seconds())
    evnts = np.array(evnts)
    intrvls = np.array(intrvls)
    print('done')

    def one(t, *args):
        return 1
    if (pmin is None) and (pmax is None):
        # поиск периодов на временах порядка часов
        hpmin = np.min(evnts[1:]-evnts[:-1])
        hpmax = 3600*24*1
        print('Searching periods from ' + str(hpmin/(3600)) +
              ' hours to ' + str(hpmax/(3600)) + ' hours')
        hper, hstat, hstat_expo = periods_statistic(
            evnts, intrvls, nbins, hpmin, hpmax, n_steps=n_steps)

        print('Simulations')
        freq = len(evnts)/np.sum(intrvls.T[1]-intrvls.T[0])
        hxstats = np.zeros((num_simulations, n_steps))
        for i in range(num_simulations):
            # генерирование пуассовонских событий
            x = Events(freq, time_intervals, one, 0)
            xevnts = np.array(x.events_in_seconds)
            xintrvls = np.array(x.intervals_in_seconds)
            print(str(i+1) + '/' + str(num_simulations))
            hxper, xstat, xstat_expo = periods_statistic(
                xevnts, intrvls, nbins, hpmin, hpmax, n_steps=n_steps)
            hxstats[i] = xstat_expo

        if saving:
            print('Saving data: ', end='done')

            hours_data = {
                'Periods': hper,
                'Chi square (без учета экспозиции)': hstat,
                'Chi square (с учетом экспозиции)': hstat_expo,
            }
            for i in range(num_simulations):
                hours_data['Simulation ' + str(i)] = hxstats[i]
            df1 = pd.DataFrame(hours_data)
            df1.to_csv(saving_directory + name + ':hours.csv', index=False)
            print('done')

        # поиск периодов на временах порядка дней
        dpmin = 3600*24*1
        dpmax = 3600*24*400
        print('Searching periods from ' + str(dpmin/(3600*24)) +
              ' day to ' + str(dpmax/(3600*24)) + ' days')
        dper, dstat, dstat_expo = periods_statistic(
            evnts, intrvls, nbins, dpmin, dpmax, n_steps=n_steps)

        print('Simulations')
        dxstats = np.zeros((num_simulations, n_steps))
        for i in range(num_simulations):
            x = Events(freq, time_intervals, one, 0)
            xevnts = np.array(x.events_in_seconds)
            xintrvls = np.array(x.intervals_in_seconds)
            print(str(i+1) + '/' + str(num_simulations))
            dxper, xstat, xstat_expo = periods_statistic(
                xevnts, intrvls, nbins, dpmin, dpmax, n_steps=n_steps)
            dxstats[i] = xstat_expo

        if saving:
            print('Saving data', end='done')
            days_data = {
                'Periods': dper,
                'Chi square (без учета экспозиции)': dstat,
                'Chi square (с учетом экспозиции)': dstat_expo,
            }
            for i in range(num_simulations):
                days_data['Simulation ' + str(i)] = dxstats[i]
            df2 = pd.DataFrame(days_data)
            df2.to_csv(saving_directory + name + ':days.csv', index=False)
            print('done')
        # поиск периодов на временах порядка года
        ypmin = 3600*24*400
        ypmax = (evnts[-1]-evnts[0])/2
        print('Searching periods from ' + str(ypmin/(3600*24*365)) +
              ' year to ' + str(ypmax/(3600*24*365)) + ' years')
        yper, ystat, ystat_expo = periods_statistic(
            evnts, intrvls, nbins, ypmin, ypmax, n_steps=n_steps)

        print('Simulations')
        yxstats = np.zeros((num_simulations, n_steps))
        for i in range(num_simulations):
            x = Events(freq, time_intervals, one, 0)
            xevnts = np.array(x.events_in_seconds)
            xintrvls = np.array(x.intervals_in_seconds)
            print(str(i+1) + '/' + str(num_simulations))
            yxper, xstat, xstat_expo = periods_statistic(
                xevnts, intrvls, nbins, ypmin, ypmax, n_steps=n_steps)
            yxstats[i] = xstat_expo

        if saving:
            print('Saving data', end='done')
            years_data = {
                'Periods': yper,
                'Chi square (без учета экспозиции)': ystat,
                'Chi square (с учетом экспозиции)': ystat_expo,
            }
            for i in range(num_simulations):
                years_data['Simulation ' + str(i)] = yxstats[i]
            df3 = pd.DataFrame(years_data)
            df3.to_csv(saving_directory + name + ':years.csv', index=False)
            print('done')
        return hper, hstat, hstat_expo, hxstats, dper, dstat, dstat_expo, dxstats, yper, ystat, ystat_expo, yxstats

    else:

        print('Searching periods from ' + str(pmin/(3600*24)) +
              ' days to ' + str(pmax/(3600*24)) + ' days')
        per, stat, stat_expo = periods_statistic(
            evnts, intrvls, nbins, pmin, pmax, n_steps=n_steps)

        print('Simulations')
        freq = len(evnts)/np.sum(intrvls.T[1]-intrvls.T[0])
        xstats = np.zeros((num_simulations, n_steps))
        for i in range(num_simulations):
            x = Events(freq, time_intervals, one, 0)
            xevnts = np.array(x.events_in_seconds)
            xintrvls = np.array(x.intervals_in_seconds)
            print(str(i+1) + '/' + str(num_simulations))
            xper, xstat, xstat_expo = periods_statistic(
                xevnts, intrvls, nbins, pmin, pmax, n_steps=n_steps)
            xstats[i] = xstat_expo

        if saving:
            print('Saving data', end='done')
            years_data = {
                'Periods': per,
                'Chi square (без учета экспозиции)': stat,
                'Chi square (с учетом экспозиции)': stat_expo,
            }
            for i in range(num_simulations):
                years_data['Simulation ' + str(i)] = xstats[i]
            df3 = pd.DataFrame(years_data)
            df3.to_csv(saving_directory + name + ':custom.csv', index=False)
            print('done')
        return per, stat, stat_expo, xstats
