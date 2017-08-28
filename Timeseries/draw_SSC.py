# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
from datetime import datetime
from pylab import * 
import math as mth
import matplotlib.dates as mdates
import matplotlib.pylab as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import FormatStrFormatter
import scipy.io as sio 

plt.rcParams['font.sans-serif']=['SimHei']
from matplotlib.pylab import rcParams

rcParams['figure.figsize'] = 10, 8

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

    theta = 29.341  # 偏转角度

    alpha2=(ts_dir2 - theta) / 180 * pi
    alpha3=(ts_dir3 - theta) / 180 * pi
    alpha4=(ts_dir4 - theta) / 180 * pi
    alpha5=(ts_dir5 - theta) / 180 * pi

    ts_vel2_ue = ts_vel2 * alpha2.map(lambda x:mth.cos(x))
    ts_vel2_un = ts_vel2 * alpha2.map(lambda x:mth.sin(x))
    
    ts_vel3_ue = ts_vel3 * alpha3.map(lambda x:mth.cos(x))
    ts_vel3_un = ts_vel3 * alpha3.map(lambda x:mth.sin(x))

    ts_vel4_ue = ts_vel4 * alpha4.map(lambda x: mth.cos(x))
    ts_vel4_un = ts_vel4 * alpha4.map(lambda x: mth.sin(x))

    ts_vel5_ue = ts_vel5 * alpha5.map(lambda x: mth.cos(x))
    ts_vel5_un = ts_vel5 * alpha5.map(lambda x: mth.sin(x))
    
    ts_flux2_ue=ts_vel2_ue*ts_ssc_25cm*0.1
    ts_flux2_un=ts_vel2_un*ts_ssc_25cm*0.1
    
    ts_ssc_avg=(ts_ssc_25cm+ts_ssc_50cm)/2
    
    ts_flux3_ue=ts_vel3_ue*ts_ssc_avg*0.1
    ts_flux3_un=ts_vel3_un*ts_ssc_avg*0.1
    
    ts_flux4_ue=ts_vel4_ue*ts_ssc_avg*0.1
    ts_flux4_un=ts_vel4_un*ts_ssc_avg*0.1
    
    ts_flux5_ue=ts_vel5_ue*ts_ssc_50cm*0.1
    ts_flux5_un=ts_vel5_un*ts_ssc_50cm*0.1    
    
    ts_allflux_ue=ts_flux2_ue+ts_flux3_ue+ts_flux4_ue+ts_flux5_ue
    ts_allflux_un=ts_flux2_un+ts_flux3_un+ts_flux4_un+ts_flux5_un

    data['flux2ue'] = ts_flux2_ue
    data['flux3ue'] = ts_flux3_ue
    data['flux4ue'] = ts_flux4_ue
    data['flux5ue'] = ts_flux5_ue

    data['flux2un'] = ts_flux2_un
    data['flux3un'] = ts_flux3_un
    data['flux4un']=ts_flux4_un
    data['flux5un']=ts_flux5_un    
    
    data['fluxue']=ts_allflux_ue
    data['fluxun'] = ts_allflux_un
                   
    writer = pd.ExcelWriter('output.xlsx')    
    data.to_excel(writer,'Sheet1')    
    writer.save()

    ts_0 = ts_surface_ssc['2014/9/25 10:00:00':'2014/9/30 6:00']
    
    font = {'family': 'serif', 'color': 'darkred', 'weight': 'normal', 'size': 6}
    
    myfont=matplotlib.font_manager.FontProperties(family='serif',size='6')  #‘serif’, ‘sans-serif’, ‘cursive’, ‘fantasy’, or ‘monospace’.

    # plt.figure(figsize=(8,10))
    fig, ax = plt.subplots(nrows=6, ncols=1, sharex=True)
    majorFormatter = FormatStrFormatter('%.2f')  # set the format of the ticker

    ax[0].plot(ts_wl, label='Depth')
    ax[0].legend(loc='upper right', fontsize='6', prop=myfont)
    ax[0].set_ylabel('Depth(m)', fontdict=font)

    ax[1].plot(ts_0, 'o-', markersize=2, label='Surface')
    ax[1].legend(loc='upper right', fontsize='2', prop=myfont)
    ax[1].set_ylabel('$SSC(kg/m^{{{3}}})$', fontdict=font)

    ax[2].plot(ts_ssc_50cm, label='50cmab')
    ax[2].legend(loc='upper right', fontsize='6',prop=myfont)
    ax[2].set_ylabel('$SSC(kg/m^{{{3}}})$', fontdict=font)

    ax[3].plot(ts_ssc_25cm, label='25cmab')
    ax[3].legend(loc='upper right', fontsize='6',prop=myfont)
    ax[3].set_ylabel('$SSC(kg/m^{{{3}}})$', fontdict=font)
    ax[3].yaxis.set_major_formatter(majorFormatter)

    ax[4].plot(ts_allflux_ue, label='along_flux')
    ax[4].legend(loc='upper right', fontsize='6',prop=myfont)
    ax[4].set_ylabel('$Flux(kg/ms)$', fontdict=font)
    ax[4].yaxis.set_major_formatter(majorFormatter)
    ax[4].axhline(y=0, color='k',linestyle='--',linewidth='0.5')

    ax[5].plot(ts_allflux_un, label='across_flux')
    ax[5].legend(loc='upper right', fontsize='6',prop=myfont)
    ax[5].set_ylabel('$Fulx(kg/ms)$', fontdict=font)
    ax[5].yaxis.set_major_formatter(majorFormatter)
    ax[5].axhline(y=0, color='k', linestyle='--', linewidth='0.5')

    ax[5].yaxis.set_major_formatter(majorFormatter)
    ax[5].xaxis.set_major_locator(days)
    ax[5].xaxis.set_major_formatter(yearsFmt)

    fig.subplots_adjust(wspace=0, hspace=0)

    for axs in ax[:]:
    #   The y-ticks will overlap with "hspace=0", so we'll hide the bottom tick
        axs.set_yticks(axs.get_yticks()[1:])
        
        for tick in axs.yaxis.get_major_ticks():
            tick.label.set_fontsize(6)
        # specify integer or one of preset strings, e.g.
        # tick.label.set_fontsize('x-small')
        # tick.label.set_rotation('vertical')
    # ax4.format_xdata = mdates.DateFormatter('%Y-%m-%d')
    fig.autofmt_xdate()
    #plt.show()
    figname = 'ts.png'    
    plt.savefig(figname, dpi = 300, bbox_inches = 'tight')
    plt.close()

if __name__ == '__main__':
    draw_ts()
