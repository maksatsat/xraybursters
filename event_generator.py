import numpy as np
import random


import numpy as np
import random


class Events:
    def __init__(self, m, intervals, period1, period2):
        self.m = m  # Число событий в единицу времени
        self.period1 = period1  # duration of events
        self.period2 = period2  # duration between events
        self.intervals = intervals
        self.start_time, self.intervals_in_seconds = Events.dates_to_seconds(
            self.intervals)
        self.allowed_periods = Events.generate_periods(
            self.period1, self.period2, self.intervals_in_seconds)
        self.events_in_seconds_initial = Events.generate_events(
            self.m, self.intervals_in_seconds, period1, period2)
        self.events_in_seconds = Events.fit_events_to_periods(
            self.events_in_seconds_initial, self.allowed_periods)
        self.events = Events.seconds_to_dates(
            self.start_time, self.events_in_seconds)

    def generate_event(m, interval, period1, period2):
        start = interval[0]
        end = interval[1]
        time = start
        events = []
        intervals = []
        while True:
            # Времена между событиями имеют экспоненциальное распределение
            x = -np.log(1.0 - random.random()) / m
            time += x
            if (time > end):
                break
            intervals.append(x)
            events.append(time)

        return events

    def generate_events(m, intervals, period1, period2):
        events = []
        for interval in intervals:
            event = Events.generate_event(m, interval, period1, period2)
            events += event
        return events

    def dates_to_seconds(intervals):
        import datetime
        start = datetime.datetime.strptime(
            intervals[0][0], '%Y-%m-%d %H:%M:%S.%f')
        intervals_in_seconds = []
        for interval in intervals:
            t_start = datetime.datetime.strptime(
                interval[0], '%Y-%m-%d %H:%M:%S.%f')-start
            t_end = datetime.datetime.strptime(
                interval[1], '%Y-%m-%d %H:%M:%S.%f')-start
            intervals_in_seconds.append(
                [t_start.total_seconds(), t_end.total_seconds()])
        return start, intervals_in_seconds

    def seconds_to_dates(start_time, events_in_seconds):
        import datetime
        event_times = []
        for event in events_in_seconds:
            event_times.append(datetime.timedelta(seconds=event)+start_time)
        return event_times

    def generate_periods(period1, period2, intervals_in_seconds):
        import random
        start, end = intervals_in_seconds[0][0], intervals_in_seconds[-1][-1]
        t = random.uniform(0, period1)
        res = []
        while t < end:
            if t+period1 < end:
                res.append([t, t+period1])
            else:
                res.append([t, end])
            t += period1 + period2
        return res

    def fit_events_to_periods(events, periods):
        res = []
        for event in events:
            for period in periods:
                if period[0] <= event <= period[1]:
                    res.append(event)
        return res


def Check_Poissoness(events):
    evts = np.array(events.events_in_seconds)
    time_intervals = evts[1:]-evts[:-1]
    time_intervals = np.delete(time_intervals, np.where(
        time_intervals > (time_intervals.std() + time_intervals.mean())))  # удаляю элементы вне 1-сигма
    m0 = events.m
    m1 = np.power(1/np.mean(time_intervals**1), 1/1)
    m2 = np.power(2/np.mean(time_intervals**2), 1/2)
    m3 = np.power(6/np.mean(time_intervals**3), 1/3)
    return m0, m1, m2, m3
