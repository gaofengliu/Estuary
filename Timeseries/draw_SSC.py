# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
from datetime import datetime
import matplotlib
import matplotlib.pylab as plt
#plt.style.use(['classic'])
plt.rcParams['font.sans-serif']=['SimHei']
from matplotlib.pylab import rcParams
rcParams['figure.figsize'] = 10, 5

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 8}

matplotlib.rc('font', **font)

matplotlib.rcParams.update({'font.size': 10})

def draw_ts():
    dateparse = lambda dates: pd.datetime.strptime(dates, '%Y/%m/%d %H:%M')
    data = pd.read_csv('test50cm_SSC.csv', parse_dates=['Time'],index_col='Time', date_parser=dateparse)
    ts_50cm_ssc=data['SSC_50cm(kg/m3)']
    ts_50cm_dep=data['depth']
    
    data = pd.read_csv('test25cm_SSC.csv', parse_dates=['Time'],index_col='Time', date_parser=dateparse)
    ts_25cm_ssc=data['SSC_25cm(kg/m3)']
    ts_25cm_dep=data['depth']
    
    ts_1=ts_50cm_ssc['2014-9-25':'2014-9-30']
    ts_2=ts_50cm_dep['2014-9-25':'2014-9-30']
    
    ts_3=ts_25cm_ssc['2014-9-25':'2014-9-30']
    ts_4=ts_25cm_dep['2014-9-25':'2014-9-30']

    plt.figure(figsize=(8,4))
    fig1, (ax1, ax2, ax3) = plt.subplots(nrows=3, ncols=1, sharex=True)
     
    #plt.subplot(411,sharex=True)
    ax1.plot(ts_1, label='SSC(50cmab)')
    ax1.legend(loc='upper right',fontsize=10)

    #plt.subplot(413,sharex=True)
    ax2.plot(ts_3,label='SSC(25cmab)')
    ax2.legend(loc='upper right',fontsize=10)

    #plt.subplot(414)
    ax3.plot(ts_4, label='Depth')
    ax3.legend(loc='upper right',fontsize=10)

    plt.subplots_adjust(wspace=0,hspace=0)
    
    for ax in [ax1, ax2]:
    #   plt.setp(ax.get_xticklabels(), visible=False)
    #   The y-ticks will overlap with "hspace=0", so we'll hide the bottom tick
        ax.set_yticks(ax.get_yticks()[1:])  
        
    

    
   # plt.tight_layout()

    plt.show()
    
    figname = 'ts.png'
    plt.savefig(figname, dpi = 300, bbox_inches = 'tight')

    plt.close()


if __name__ == '__main__':
    draw_ts()
