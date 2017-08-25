# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
from datetime import datetime
import matplotlib.dates as mdates
import matplotlib.pylab as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import FormatStrFormatter
import scipy.io as sio 

plt.rcParams['font.sans-serif']=['SimHei']
from matplotlib.pylab import rcParams
rcParams['figure.figsize'] = 10, 5

def draw_ts():
    years = mdates.YearLocator()  # every year
    months = mdates.MonthLocator()  # every month
    days = mdates.DayLocator()  # every day
    yearsFmt = mdates.DateFormatter('%Y-%m-%d')

    dateparse = lambda dates: pd.datetime.strptime(dates, '%Y/%m/%d %H:%M')
    data = pd.read_csv('20141008_YJG.csv', parse_dates=['DateTime'], index_col='DateTime', date_parser=dateparse)

    data2 = pd.read_csv('surface_SSC.csv', parse_dates=['Time'], index_col='Time', date_parser=dateparse)
    ts_surface_ssc = data2['SSC_surface(kg/m3)']

    ts_wl = data['Pressure']
    ts_vel2 = data['Speed#1(0.2m)']
    ts_dir2 = data['Dir#1(0.2m)']
    ts_vel3 = data['Speed#2(0.3m)']
    ts_dir3 = data['Dir#2(0.3m)']
    ts_vel4 = data['Speed#3(0.4m)']
    ts_dir4 = data['Dir#3(0.4m)']
    ts_vel5 = data['Speed#4(0.5m)']
    ts_dir5 = data['Dir#4(0.5m)']

    ts_ssc_50cm = data['SSC_50cm']
    ts_ssc_25cm = data['SSC_25cm']

    pi = 3.1415926

    theta = 29.341 / 360 * pi

    # matlab文件名
    # matfn=u'F:/0python/0python/0professional code/Estuary/Timeseries/[20141008]AquaPro_OBS3+_CH1_T8913_25cmab_CH2_T8922_50cmab'
    # data=sio.loadmat(matfn)

    ts_0 = ts_surface_ssc['2014/9/25 10:00:00':'2014/9/30 6:00']
    font = {'family': 'serif', 'color': 'darkred', 'weight': 'normal', 'size': 8}

    plt.figure(figsize=(8,4))
    fig, ax = plt.subplots(nrows=4, ncols=1, sharex=True)

    majorFormatter = FormatStrFormatter('%.2f')  # set the format of the ticker

    ax[0].plot(ts_0, 'o-', markersize=5, label='Surface')
    ax[0].legend(loc='upper right', fontsize='10')
    # ax1.tick_params(direction='in', length=6, width=2, colors='k',labelleft='on')
    ax[0].set_ylabel('$SSC(kg/m^{{{3}}})$', fontdict=font)

    ax[1].plot(ts_ssc_50cm, label='50cmab')
    ax[1].legend(loc='upper right', fontsize='10')
    ax[1].set_ylabel('$SSC(kg/m^{{{3}}})$', fontdict=font)

    ax[2].plot(ts_ssc_25cm, label='25cmab')
    ax[2].legend(loc='upper right', fontsize='10')
    ax[2].set_ylabel('$SSC(kg/m^{{{3}}})$', fontdict=font)
    ax[2].yaxis.set_major_formatter(majorFormatter)

    ax[3].plot(ts_wl, label='Depth')
    ax[3].legend(loc='upper right', fontsize='10')
    ax[3].set_ylabel('Depth(m)', fontdict=font)
    ax[3].yaxis.set_major_formatter(majorFormatter)

    ax[3].xaxis.set_major_locator(days)
    ax[3].xaxis.set_major_formatter(yearsFmt)

    fig.subplots_adjust(wspace=0, hspace=0.15)

    for axs in ax[:]:
    #   The y-ticks will overlap with "hspace=0", so we'll hide the bottom tick
        axs.set_yticks(axs.get_yticks()[1:])
    for tick in axs.yaxis.get_major_ticks():
        tick.label.set_fontsize(8)
        # specify integer or one of preset strings, e.g.
        # tick.label.set_fontsize('x-small')
        # tick.label.set_rotation('vertical')
    # ax4.format_xdata = mdates.DateFormatter('%Y-%m-%d')
    fig.autofmt_xdate()
    plt.show()
    figname = 'ts.png'
    plt.savefig(figname, dpi = 300, bbox_inches = 'tight')
    plt.close()



if __name__ == '__main__':
    draw_ts()
