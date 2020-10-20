import math
import random


class Events:
    def __init__(self, m, intervals):
        self.m = m  # Число событий в единицу времени
        self.intervals = intervals
        self.events = Events.generate_events(self.m, self.intervals)

    def generate_event(m, interval):
        start = interval[0]
        end = interval[1]
        time = start
        events = []
        intervals = []
        while True:
            # Времена между событиями имеют экспоненциальное распределение
            x = -math.log(1.0 - random.random()) / m
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
