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
years = mdates.YearLocator()  # every year
months = mdates.MonthLocator()  # every month
days = mdates.DayLocator()  # every day
yearsFmt = mdates.DateFormatter('%Y-%m-%d')

dateparse = lambda dates: pd.datetime.strptime(dates, '%Y/%m/%d %H:%M')
data = pd.read_csv('20141008_YJG.csv', parse_dates=['DateTime'], index_col='DateTime', date_parser=dateparse)
data2 = pd.read_csv('surface_SSC.csv', parse_dates=['Time'], index_col='Time', date_parser=dateparse)

ts_surface_ssc = data2['SSC_surface(kg/m3)']
ts_wl = data['Pressure']
col=data.columns
ts_vel=data[col[1:9:2]]
ts_dir=data[col[2:10:2]]
ts_ssc_50cm = data['SSC_50cm']
ts_ssc_25cm = data['SSC_25cm']

pi = 3.1415926
theta = 29.341 
alpha=(ts_dir - theta) / 180 * pi
alpha.columns =[colname for colname in ts_vel.columns]

ts_vel_ue=ts_vel*alpha.applymap(lambda x:mth.cos(x))
ts_vel_un=ts_vel*alpha.applymap(lambda x:mth.sin(x))

ts_ssc_avg=(ts_ssc_25cm+ts_ssc_50cm)/2
ts_ssc=ts_vel
velcol=ts_vel.columns

ts_ssc.loc[:,velcol[0]]=ts_ssc_25cm
ts_ssc.loc[:,velcol[1]]=ts_ssc_avg
ts_ssc.loc[:,velcol[2]]=ts_ssc_avg
ts_ssc.loc[:,velcol[3]]=ts_ssc_50cm

ts_flux_ue=ts_vel_ue*ts_ssc*0.1
ts_flux_un=ts_vel_un*ts_ssc*0.1

ts_flux_ue.columns = ['flux2ue', 'flux3ue', 'flux4ue', 'flux5ue']
ts_flux_un.columns = ['flux2un', 'flux3un', 'flux4un', 'flux5un']

ts_flux_ue['along_flux']=ts_flux_ue.sum(1)   # 存在的问题是 nan相加得到了0
ts_flux_un['across_flux']=ts_flux_un.sum(1)

writer = pd.ExcelWriter('output2.xlsx')
ts_flux_ue.to_excel(writer,'Sheet1')
ts_flux_un.to_excel(writer,'Sheet2')
writer.save()

ts_0 = ts_surface_ssc['2014/9/25 10:00:00':'2014/9/30 6:00']    
font = {'family': 'serif', 'color': 'darkred', 'weight': 'normal', 'size': 6}
myfont=matplotlib.font_manager.FontProperties(family='serif',size='6')  #‘serif’, ‘sans-serif’, ‘cursive’, ‘fantasy’, or ‘monospace’.

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

ax[4].plot(ts_flux_ue['along_flux'], label='along_flux')
ax[4].legend(loc='upper right', fontsize='6',prop=myfont)
ax[4].set_ylabel('$Flux(kg/ms)$', fontdict=font)
ax[4].yaxis.set_major_formatter(majorFormatter)
ax[4].axhline(y=0, color='k',linestyle='--',linewidth='0.5')

ax[5].plot(ts_flux_un['across_flux'], label='across_flux')
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
figname = 'ts2.png'    
plt.savefig(figname, dpi = 300, bbox_inches = 'tight')
plt.close()
