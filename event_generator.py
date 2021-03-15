import numpy as np
import random


class Events:
    def __init__(self, m, intervals, shape_function, period):
        self.m = m  # Число событий в единицу времени
        self.period = period
        self.shape_function = shape_function
        self.intervals = intervals
        self.start_time, self.intervals_in_seconds = Events.dates_to_seconds(
            self.intervals)
        self.events_in_seconds_initial = Events.generate_events(
            self.m, self.intervals_in_seconds, self.period)
        self.events_in_seconds = Events.signal_shape(
            self.events_in_seconds_initial, self.shape_function, self.period)
        self.events = Events.seconds_to_dates(
            self.start_time, self.events_in_seconds)

    def generate_events_for_interval(m, interval, period):
        time = interval[0] - np.log(1.0 - random.random()) / m
        events = []
        while time < interval[1]:
            events.append(time)
            time += -np.log(1.0 - random.random()) / m
            # Времена между событиями имеют экспоненциальное распределение
        return events

    def generate_events(m, intervals, period):
        events = []
        for interval in intervals:
            event = Events.generate_events_for_interval(m, interval, period)
            events += event  # ?
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

    # def generate_periods(period, intervals_in_seconds):
    #     import random
    #     start, end = intervals_in_seconds[0][0], intervals_in_seconds[-1][-1]
    #     t = random.uniform(0, period1)
    #     res = []
    #     while t < end:
    #         if t+period1 < end:
    #             res.append([t, t+period1])
    #         else:
    #             res.append([t, end])
    #         t += period1 + period2
    #     return res

    def fit_events_to_periods(events, periods):
        res = []
        for event in events:
            for period in periods:
                if period[0] <= event <= period[1]:
                    res.append(event)
        return res

    def signal_shape(events, shape_function, period):
        res = []
        for event in events:
            if np.random.random() < shape_function(event, period):
                res.append(event)
        return res

    def sin_signal(t, period, phase=0, m=0.5):
        # period = period1 + period2
        return 1-m*(1-np.sin(2*np.pi*t/period+phase))

    def pulse_wave(t, period1, period2, phase=0):
        period1 = period
        period2 = 0.5 * period
        if (t-phase) % (period1+period2) > period1:
            return 1
        else:
            return 0


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
