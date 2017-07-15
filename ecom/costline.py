# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 18:53:36 2016

@author: admin
"""
import numpy as np
import struct
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os
import datetime as dt
from scipy.interpolate import interp1d
from scipy.interpolate import griddata


indir = 'LandXY/'

findir = os.listdir(indir)
island = [x for x in range(len(findir))]
 
for i in range(len(findir)):
    island[i]=np.loadtxt(indir+findir[i], usecols=(0,1) , dtype=float) # float64 narray


plt.figure(figsize=(10,10),dpi=600)

ax=plt.subplot(221)

for i in range(len(findir)):    
    plt.plot(island[i][:,0],island[i][:,1], color="black", linewidth=0.8, linestyle="-")



xmin=306093.*0.9
xmax=469822.*1.1
ymin=3399766.*0.98
ymax=3539519.*1.02

xlm=np.linspace(xmin,xmax,5,endpoint=True)
xtik=xlm/1000.

plt.xlim(xmin,xmax)
plt.xticks(xlm,xtik)

for label in ax.xaxis.get_ticklabels():
    # label is a Text instance
    label.set_color('red')
    label.set_rotation(45)
    label.set_fontsize(8)

formatter = ticker.FormatStrFormatter('%1.1f')
ax.xaxis.set_major_formatter(formatter)

ylm=np.linspace(ymin,ymax,5,endpoint=True)
ytik=ylm/1000.



plt.ylim(ymin,ymax)
plt.yticks(ylm,ytik)

for label in ax.yaxis.get_ticklabels():
    # label is a Text instance
    label.set_color('red')
    label.set_fontsize(8)

ax.yaxis.set_major_formatter(formatter)


plt.xlabel('X(km)')
plt.ylabel('Y(km)')
# plt.style.use('bmh')

figname = 'test'+ '.png'
plt.savefig(figname, dpi = 600, bbox_inches = 'tight')


plt.close()