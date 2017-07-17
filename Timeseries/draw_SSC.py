# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
from datetime import datetime
import matplotlib.dates as mdates
import matplotlib.pylab as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import FormatStrFormatter

plt.rcParams['font.sans-serif']=['SimHei']
from matplotlib.pylab import rcParams
rcParams['figure.figsize'] = 10, 5



def draw_ts():
    years = mdates.YearLocator()  # every year
    months = mdates.MonthLocator()  # every month
    days = mdates.DayLocator()  # every day
    yearsFmt = mdates.DateFormatter('%Y-%m-%d')

    dateparse = lambda dates: pd.datetime.strptime(dates, '%Y/%m/%d %H:%M')
    data = pd.read_csv('test50cm_SSC.csv', parse_dates=['Time'],index_col='Time', date_parser=dateparse)
    ts_50cm_ssc=data['SSC_50cm(kg/m3)']
    ts_50cm_dep=data['depth']

    data = pd.read_csv('test25cm_SSC.csv', parse_dates=['Time'],index_col='Time', date_parser=dateparse)
    ts_25cm_ssc=data['SSC_25cm(kg/m3)']
    ts_25cm_dep = data['depth2']

    data = pd.read_csv('test_surface_SSC.csv', parse_dates=['Time'], index_col='Time', date_parser=dateparse)
    ts_surface_ssc = data['SSC_surface(kg/m3)']

    ts_1 = ts_50cm_ssc['2014-9-25':'2014-10-3']
    ts_2 = ts_50cm_dep['2014-9-25':'2014-10-3']

    ts_3 = ts_25cm_ssc['2014-9-25':'2014-10-3']
    ts_4 = ts_25cm_dep['2014-9-25':'2014-10-3']

    ts_0 = ts_surface_ssc['2014/9/25  10:00:00':'2014/9/30 6:00']

    font = {'family': 'serif',
            'color': 'darkred',
            'weight': 'normal',
            'size': 8}

    plt.figure(figsize=(8,4))

    fig, ax = plt.subplots(nrows=4, ncols=1, sharex=True)

    majorFormatter = FormatStrFormatter('%.2f')  # set the format of the ticker

    ax[0].plot(ts_0, 'o-', markersize=5, label='Surface')
    ax[0].legend(loc='upper right', fontsize='10')
    # ax1.tick_params(direction='in', length=6, width=2, colors='k',labelleft='on')
    ax[0].set_ylabel('$SSC(kg/m^{{{3}}})$', fontdict=font)

    ax[1].plot(ts_1, label='50cmab')
    ax[1].legend(loc='upper right', fontsize='10')
    ax[1].set_ylabel('$SSC(kg/m^{{{3}}})$', fontdict=font)

    ax[2].plot(ts_3, label='25cmab')
    ax[2].legend(loc='upper right', fontsize='10')
    ax[2].set_ylabel('$SSC(kg/m^{{{3}}})$', fontdict=font)
    ax[2].yaxis.set_major_formatter(majorFormatter)

    ax[3].plot(ts_4, label='Depth')
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
