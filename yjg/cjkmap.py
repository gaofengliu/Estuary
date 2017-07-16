__author__ = 'mac'
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pylab as plt
from matplotlib.pylab import rcParams

rcParams['figure.figsize']= 10,4
font = {'family': 'normal','weight':'bold','size':10}

plt.rcParams['font.sans-serif']=['SimHei']

matplotlib.rc('font',**font)

def draw_cjk_map():
#    with open('yjg_coastline.txt') as f:
#         coastlines = f.readlines()
    yjgcoast = np.loadtxt('yjg_coastline.txt')
    hbykmt = np.loadtxt('hbykmt.txt')
    hjmt = np.loadtxt('hjmt.txt')
    lycmt = np.loadtxt('lycmt.txt')
    shmt = np.loadtxt('shmt.txt')

    fig,ax=plt.subplots()    
    
    ax.plot(yjgcoast[:,0],yjgcoast[:,1],color='k',linewidth='2') 
    ax.plot(hbykmt[:,0],hbykmt[:,1],color='k')
    ax.plot(hjmt[:,0],hjmt[:,1],color='k') 
    ax.plot(lycmt[:,0],lycmt[:,1],color='k')
    ax.plot(shmt[:,0],shmt[:,1],color='k')
    
    ax.get_yaxis().get_major_formatter().set_useOffset(False)

    


#    startline = 0
#    for ii in range(startline,startline+175):
#        data[ii] = coastlines[ii].split()







if __name__ == '__main__':
    draw_cjk_map()


