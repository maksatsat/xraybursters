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


def efsearch(bursts_url, osbervations_url, bins=10, n_steps=1000, num_simulations=100):

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
    name = bursts_table['name'][0]
    instruments = bursts_table['instr'].unique()
    print('done')
    print('Mergin observation intervals: ', end='')
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
    print('Merging burst times', end='')
    # объединенные времена вспышек

    def notfoundsim(b, merged_bursts):
        # if len(merged_bursts)==0:
            # return True
        for x in merged_bursts:
            delta = (Time(x, format='mjd')-Time(b, format='mjd')
                     ).to_datetime().total_seconds()
            if (np.abs(delta) <= 30):
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
    print('done')
    time_intervals = []
    for interv in merged_intervals:
        tm1 = Time(interv[0], format='mjd')
        d1 = datetime.datetime.strptime(
            tm1.iso, '%Y-%m-%d %H:%M:%S.%f').strftime("%Y-%m-%d %H:%M:%S.%f")
        tm2 = Time(interv[1], format='mjd')
        d2 = datetime.datetime.strptime(
            tm2.iso, '%Y-%m-%d %H:%M:%S.%f').strftime("%Y-%m-%d %H:%M:%S.%f")
        time_intervals.append([d1, d2])

    events = []
    for burst in merged_bursts:
        tm = Time(burst, format='mjd')
        d = datetime.datetime.strptime(tm.iso, '%Y-%m-%d %H:%M:%S.%f')
        events.append(d)

    start, intrvls = Events.dates_to_seconds(time_intervals)
    evnts = []
    for event in events:
        evnts.append((event-start).total_seconds())
    evnts = np.array(evnts)
    intrvls = np.array(intrvls)
    print('done')

    # поиск периодов на временах порядка часов
    pmin = (events[1]-events[0]).total_seconds()
    pmax = 3600*24*1
    print('Searching periods from ' + str(pmin/(3600)) +
          ' hours to ' + str(pmax/(3600)) + ' hours')
    hper, hstat, hstat_expo = periods_statistic(
        evnts, intrvls, bins, pmin, pmax, n_steps=n_steps)

    def one(t, *args):
        return 1

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
            xevnts, intrvls, bins, pmin, pmax, n_steps=n_steps)
        hxstats[i] = xstat_expo

    # поиск периодов на временах порядка дней
    pmin = 3600*24*1
    pmax = 3600*24*400
    print('Searching periods from ' + str(pmin/(3600*24)) +
          ' days to ' + str(pmax/(3600*24)) + ' days')
    dper, dstat, dstat_expo = periods_statistic(
        evnts, intrvls, bins, pmin, pmax, n_steps=n_steps)

    print('Simulations')
    dxstats = np.zeros((num_simulations, n_steps))
    for i in range(num_simulations):
        x = Events(freq, time_intervals, one, 0)
        xevnts = np.array(x.events_in_seconds)
        xintrvls = np.array(x.intervals_in_seconds)
        print(str(i+1) + '/' + str(num_simulations))
        dxper, xstat, xstat_expo = periods_statistic(
            xevnts, intrvls, bins, pmin, pmax, n_steps=n_steps)
        dxstats[i] = xstat_expo

    # поиск периодов на временах порядка года
    pmin = 3600*24*400
    pmax = (evnts[-1]-evnts[0])/2
    print('Searching periods from ' + str(pmin/(3600*24*365)) +
          ' days to ' + str(pmax/(3600*24*365)) + ' days')
    yper, ystat, ystat_expo = periods_statistic(
        evnts, intrvls, bins, pmin, pmax, n_steps=n_steps)

    print('Simulations')
    yxstats = np.zeros((num_simulations, n_steps))
    for i in range(num_simulations):
        x = Events(freq, time_intervals, one, 0)
        xevnts = np.array(x.events_in_seconds)
        xintrvls = np.array(x.intervals_in_seconds)
        print(str(i+1) + '/' + str(num_simulations))
        yxper, xstat, xstat_expo = periods_statistic(
            xevnts, intrvls, bins, pmin, pmax, n_steps=n_steps)
        yxstats[i] = xstat_expo

    print('Saving data')

    hours_data = {
        'Periods': hper,
        'Chi square (без учета экспозиции)': hstat,
        'Chi square (с учетом экспозиции)': hstat,
    }
    for i in range(num_simulations):
        hours_data['Simulation ' + str(i)] = hxstats[i]
    df1 = pd.DataFrame(hours_data)
    df1.to_csv('results/' + name + '/' + name + ' :hours.csv', index=False)

    days_data = {
        'Periods': dper,
        'Chi square (без учета экспозиции)': dstat,
        'Chi square (с учетом экспозиции)': dstat,
    }
    for i in range(num_simulations):
        days_data['Simulation ' + str(i)] = dxstats[i]
    df2 = pd.DataFrame(days_data)
    df2.to_csv('results/' + name + '/' + name + ' :days.csv', index=False)

    years_data = {
        'Periods': yper,
        'Chi square (без учета экспозиции)': ystat,
        'Chi square (с учетом экспозиции)': ystat,
    }
    for i in range(num_simulations):
        years_data['Simulation ' + str(i)] = yxstats[i]
    df3 = pd.DataFrame(years_data)
    df3.to_csv('results/' + name + '/' + name + ' :years.csv', index=False)

    return hper, hstat, hstat_expo, hxstats, dper, dstat, dstat_expo, dxstats, yper, ystat, ystat_expo, yxstats
