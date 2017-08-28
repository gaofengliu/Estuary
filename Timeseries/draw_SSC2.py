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
rcParams['figure.figsize'] = 10, 5
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

writer = pd.ExcelWriter('output2.xlsx')    
ts_flux_ue.to_excel(writer,'Sheet1')
ts_flux_un.to_excel(writer, 'Sheet2')

# ts_flux= pd.concat(ts_flux_ue,ts_flux_un)
writer.save()


 

       
               
