# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import math as mth
from datetime import datetime
import matplotlib.dates as mdates
import matplotlib.pylab as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import FormatStrFormatter
import scipy.io as sio

plt.rcParams['font.sans-serif'] = ['SimHei']
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
theta = 29.341

alpha2=(ts_dir2 - theta) / 180 * pi
alpha3=(ts_dir3 - theta) / 180 * pi
alpha4=(ts_dir4 - theta) / 180 * pi
alpha5=(ts_dir5 - theta) / 180 * pi

ts_vel2_ue = ts_vel2 * alpha2.map(lambda x:mth.cos(x))
ts_vel2_un = ts_vel2 * alpha2.map(lambda x:mth.sin(x))

ts_vel3_ue = ts_vel3 * alpha3.map(lambda x:mth.cos(x))
ts_vel3_un = ts_vel3 * alpha3.map(lambda x:mth.sin(x))

ts_vel4_ue = ts_vel2 * alpha4.map(lambda x:mth.cos(x))
ts_vel4_un = ts_vel2 * alpha4.map(lambda x:mth.sin(x))

ts_vel5_ue = ts_vel2 * alpha5.map(lambda x:mth.cos(x))
ts_vel5_un = ts_vel2 * alpha5.map(lambda x:mth.sin(x))

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

data['flux2ue']=ts_flux2_ue
data['flux2un']=ts_flux2_un

data['flux3ue']=ts_flux3_ue
data['flux3un']=ts_flux3_un

data['flux4ue']=ts_flux4_ue
data['flux4un']=ts_flux4_un

data['flux5ue']=ts_flux5_ue
data['flux5un']=ts_flux5_un


data['fluxue']=ts_allflux_ue
data['fluxun']=ts_allflux_un
     
     
               
writer = pd.ExcelWriter('output.xlsx')

data.to_excel(writer,'Sheet1')

writer.save()



