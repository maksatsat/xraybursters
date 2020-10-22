import numpy as np
import random


class Events:
    def __init__(self, m, intervals):
        self.m = m  # Число событий в единицу времени
        self.intervals = intervals
        self.start_time, self.intervals_in_seconds = Events.dates_to_seconds(
            self.intervals)
        self.events_in_seconds = Events.generate_events(
            self.m, self.intervals_in_seconds)
        self.events = Events.seconds_to_dates(
            self.start_time, self.events_in_seconds)

    def generate_event(m, interval):
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

    def generate_events(m, intervals):
        events = []
        for interval in intervals:
            event = Events.generate_event(m, interval)
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
