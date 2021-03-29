import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

name = input('Print name of the system\n')

flag = None
print('Choose time interval to plot: 1 - hours, 2 - days, 3 - years')
while flag is None:
    flag = input()
    if flag == '1':
        filename = name + ' :hours.csv'
        time_norm = 3600
        xlabel = 'days'
    elif flag == '2':
        filename = name + ' :days.csv'
        time_norm = 3600*24
    elif flag == '3':
        filename = name + ' :years.csv'
        time_norm = 3600*24*365
    else:
        flag = None

df = pd.read_csv(filename)

periods = df['Periods']
stats = df['Chi square (без учета экспозиции)']
stats_expo = df['Chi square (с учетом экспозиции)']
plt.plot(periods/time_norm, stats, 'r')
plt.plot(periods/time_norm, stats_expo, 'b')

plt.yscale('log')

plt.show()
