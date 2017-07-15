# -*- coding: utf-8 -*-
import numpy as np
import struct
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from scipy.interpolate import griddata
# from mpl_toolkits.axes_grid1 import make_axes_locatable
import gsw
import os
import datetime as dt
from scipy.interpolate import interp1d

case = 'case2011'
indir = '../../' + case + '/Ecomsi3dsed/output/field_distri/'
kind1 = "temperature/t_field_"
kind2 = 'salinity/s_field_'
kind3 = 'current/v_field_'
kind4 = 'elevation/el_field_'
kind5 = 'momentum/mmt_field_'


ntime1 = 744
ntime2 = ntime1+1
ntime = ntime2-ntime1+1
factor = 1.0e5
date0 = dt.datetime(2011,7,1)
grav = 9.801

date = date0 + dt.timedelta(hours=ntime1)
tracer = 0
total_tracers = 2
pressure = 0.
lon_1 = 118.05
lon_2 = 122.75
lat_1 = 28.05
lat_2 = 32.45
stick_interval_x = 1
stick_interval_y = 1

ubig = 0.5
deplevels = np.array([10.,20.,30.,40.,50.,60.,70.,80.,90.,100.,120.,150.,200.,300.,500])
[lon,lat]=np.meshgrid(np.arange(lon_1,lon_2,0.25),np.arange(lat_1,lat_2,0.25))
if not os.path.exists(case):
    os.makedirs(case)

def slice(X,Y,H,ZZ,S,zlevel):
    JM,IM = H.shape
    sz = np.zeros([JM,IM])
    if zlevel == 0:
        sz=S[0,:]
    elif zlevel == -1:
        sz=S[-2,:]
    else:
        zlevel=-1.*zlevel
        for j in range(0,JM):
            for i in range(0,IM):
                dep = H[j,i] * ZZ
                if zlevel > dep[-1] :
                    sz[j,i] = interp1d(dep,S[:-1,j,i],kind='cubic')(zlevel)
    return sz


def sphere_dis(lon1,lat1,lon2,lat2):
    pi = 3.1415926 / 180.
    rearth = 6370000.
    dx = (lon2 - lon1) * pi
    ymean = 0.5 * (lat1 + lat2) * pi
    dx = rearth * dx * np.cos(ymean)
    dy = rearth * (lat2 - lat1) * pi
    dis = np.sqrt(dx **2 + dy **2)
    return dis

#  read grid file
with open('../../grid_info/ch_hzbc.grd') as f:
    lines = f.readlines()

KB = int(lines[2].split()[0])
KBM1 = KB-1
Z = np.array([float(lines[3+k].split()[0]) for k in range(0,KB)]).squeeze()
IM = int(lines[4+KB].split()[0])
JM = int(lines[4+KB].split()[1])
ZZ = 0.5 * (Z[:-1] + Z[1:])
DZ = Z[:-1] - Z[1:]

IMJM = IM*JM
x = np.zeros([JM,IM])
y = np.zeros([JM,IM])
h = np.zeros([JM,IM])


t = np.zeros([ntime,KB,JM,IM])
s = np.zeros([ntime,KB,JM,IM])
u = np.zeros([ntime,KB,JM,IM])
v = np.zeros([ntime,KB,JM,IM])
el = np.zeros([ntime,JM,IM])


startline = 5+KB
for ii in range(startline,startline+IMJM):
    data = lines[ii].split()
    i = int(data[0])-1
    j = int(data[1])-1
    x[j,i] = data[8]
    y[j,i] = data[9]

cori_3d = np.tile(y[np.newaxis,:,:], (KB,1,1))
cori_3d = 2. * 2. * 3.1415926 / (86400.) * np.sin(cori_3d * 3.1415926 / 180.)

# read grid depth
with open('../../grid_info/ch_hzbc_dep.dat') as f:
    deplines = f.readlines()

startline = 1
for ii in range(startline,startline+IMJM):
    data = deplines[ii].split()
    i = int(data[0])-1
    j = int(data[1])-1
    h[j,i] = data[2]    

## Reading the data field
for nt in range(ntime1,ntime2+1):
    print nt
    ii = nt-ntime1
    tfile = indir + kind1 + np.str(nt).zfill(6)
    sfile = indir + kind2 + np.str(nt).zfill(6)    
    vfile = indir + kind3 + np.str(nt).zfill(6)
    efile = indir + kind4 + np.str(nt).zfill(6)
    with open (sfile,'rb') as f:
        f.read(4)
        thour = f.read(4)
        thour = np.array(struct.unpack('>f',thour))
        data = f.read(IM*JM*KB*4)
        s[ii,:]  = np.array(struct.unpack('>'+np.str(IM*JM*KB) + 'f', data)).reshape(KB,JM,IM)
        f.read(4)

    with open (vfile,'rb') as f:
        f.read(4)
        thour = f.read(4)
        thour = np.array(struct.unpack('>f',thour))
        data = f.read(IM*JM*KB*4)
        u[ii,:]  = np.array(struct.unpack('>'+np.str(IM*JM*KB) + 'f', data)).reshape(KB,JM,IM)
        data = f.read(IM*JM*KB*4)
        v[ii,:]  = np.array(struct.unpack('>'+np.str(IM*JM*KB) + 'f', data)).reshape(KB,JM,IM)
        f.read(4)
        
    with open (efile,'rb') as f:
        f.read(4)
        thour = f.read(4)
        thour = np.array(struct.unpack('>f',thour))
        data = f.read(IM*JM*4)
        el[ii,:]  = np.array(struct.unpack('>'+np.str(IM*JM) + 'f', data)).reshape(JM,IM)
        
        
tmean = 20
smean = s.mean(axis=0)
umean = u.mean(axis=0)
vmean = v.mean(axis=0)
emean = el.mean(axis=0)
rhomean = gsw.rho(smean,tmean,pressure) 
h3d = np.tile(h[np.newaxis,:],[KB,1,1])
rhomean = np.ma.masked_array(rhomean,h3d<3.)
rho0 = rhomean.mean()

plt.figure(figsize=(10,10))
print 'subplot-1 '
plt.subplot(221)
cmap = plt.cm.gist_ncar
data = emean
#vmin = data[II].min()
#vmax = data[II].max()
vmin = 0.01
vmax = 5.0
data[data<=vmin]=vmin
data[data>=vmax]=vmax
levels = np.arange(vmin,vmax+(vmax-vmin)/30.,(vmax-vmin)/30.)
m = Basemap(llcrnrlon=lon_1,llcrnrlat=lat_1,urcrnrlon=lon_2,urcrnrlat=lat_2,\
                    resolution='h',projection='merc')
m.drawparallels(np.arange(lat_1,lat_2,stick_interval_x),labels=[1,0,0,0],linewidth=0.0)
m.drawmeridians(np.arange(lon_1,lon_2,stick_interval_y),labels=[0,0,0,1],linewidth=0.0)
plt.tick_params(axis='both',labelsize=8)
m.contourf(x,y,data,levels,cmap = cmap)
#m.colorbar()
m.contour(x,y,h,deplevels,linewidths = 0.2)

plt_zlevel = 3.
um = umean[0,:]
vm = vmean[0,:]
u_plt = np.ma.masked_array(um,mask = (h*(-1.*ZZ[-1])<plt_zlevel))
v_plt = np.ma.masked_array(vm,mask = (h*(-1.*ZZ[-1])<plt_zlevel))
xx = x[~np.isnan(u_plt)]
yy = y[~np.isnan(u_plt)]
uu = u_plt[~np.isnan(u_plt)]
vv = v_plt[~np.isnan(u_plt)]
u_plt = griddata((xx,yy),uu,(lon,lat),method='cubic')
v_plt = griddata((xx,yy),vv,(lon,lat),method='cubic')
uabs = np.sqrt(u_plt**2 + v_plt**2)
u_plt[uabs<0.001] = np.nan
v_plt[uabs<0.001] = np.nan
uv = np.sqrt(u_plt*u_plt + v_plt*v_plt)
Q = m.quiver(lon[uv>=ubig],lat[uv>=ubig],u_plt[uv>=ubig],v_plt[uv>=ubig],scale=8,color='purple')
qk = plt.quiverkey(Q,0.2,0.5,0.5,'0.5m/s',color='purple') #,fontsize = 6)
Q = m.quiver(lon[uv<ubig],lat[uv<ubig],u_plt[uv<ubig],v_plt[uv<ubig],scale=3,color='black')
qk = plt.quiverkey(Q,0.2,0.4,0.2,'0.2m/s',color='black') #,fontsize = 6)

m.drawcoastlines(linewidth = 0.2)          
# m.fillcontinents(color='gray')
plt.annotate('A',xy=(0.15,0.6),xycoords='axes fraction') #,fontsize=14)

'''
print 'subplot 2'
plt.subplot(222)

z_level = 20.
height = np.zeros([JM,IM])
for i in range(0,IM):
    for j in range(0,JM):
        if h[j,i] >= -99999.:
            z_base = min(z_level,h[j,i])
            for k in range(1,KB):
                depth_k = -1. * ((h[j,i] + emean[j,i])*Z[k] + emean[j,i])
                depth_km1 = -1. * ((h[j,i] + emean[j,i])*Z[k-1] + emean[j,i])
                if depth_k <= z_base:
                    height[j,i] = height[j,i] + rhomean[k-1,j,i] * DZ[k-1] * (h[j,i] + emean[j,i]) * grav
                elif (depth_k > z_base) & (depth_km1 <= z_base):
                    ratio = (z_base - depth_km1) / (depth_k - depth_km1)
                    height[j,i] = height[j,i] + rhomean[k-1,j,i] * DZ[k-1] * (h[j,i] + emean[j,i]) * ratio * grav
            height[j,i] = height[j,i] /10000.
data = height
vmin = 0.01
vmax = 5.0
data[data<=vmin]=vmin
data[data>=vmax]=vmax
data = np.ma.masked_array(data,h<z_level)
levels = np.arange(vmin,vmax+(vmax-vmin)/30.,(vmax-vmin)/30.)
m = Basemap(llcrnrlon=lon_1,llcrnrlat=lat_1,urcrnrlon=lon_2,urcrnrlat=lat_2,\
                    resolution='h',projection='merc')
m.drawparallels(np.arange(lat_1,lat_2,stick_interval_x),labels=[1,0,0,0],linewidth=0.0)
m.drawmeridians(np.arange(lon_1,lon_2,stick_interval_y),labels=[0,0,0,1],linewidth=0.0)
plt.tick_params(axis='both',labelsize=8)
m.contourf(x,y,data,levels,cmap = cmap)
#m.colorbar()
m.contour(x,y,h,deplevels,linewidths = 0.2)
plt_zlevel = z_level
um = slice(x,y,h,ZZ,umean,plt_zlevel)
vm = slice(x,y,h,ZZ,vmean,plt_zlevel)
u_plt = np.ma.masked_array(um,mask = (h*(-1.*ZZ[-1])<plt_zlevel))
v_plt = np.ma.masked_array(vm,mask = (h*(-1.*ZZ[-1])<plt_zlevel))
xx = x[~np.isnan(u_plt)]
yy = y[~np.isnan(u_plt)]
uu = u_plt[~np.isnan(u_plt)]
vv = v_plt[~np.isnan(u_plt)]
u_plt = griddata((xx,yy),uu,(lon,lat),method='cubic')
v_plt = griddata((xx,yy),vv,(lon,lat),method='cubic')
uabs = np.sqrt(u_plt**2 + v_plt**2)
u_plt[uabs<0.001] = np.nan
v_plt[uabs<0.001] = np.nan
uv = np.sqrt(u_plt*u_plt + v_plt*v_plt)
Q = m.quiver(lon[uv>=ubig],lat[uv>=ubig],u_plt[uv>=ubig],v_plt[uv>=ubig],scale=8,color='purple')
qk = plt.quiverkey(Q,0.2,0.5,0.5,'0.5m/s',color='purple') #,fontsize = 6)
Q = m.quiver(lon[uv<ubig],lat[uv<ubig],u_plt[uv<ubig],v_plt[uv<ubig],scale=3,color='black')
qk = plt.quiverkey(Q,0.2,0.4,0.2,'0.2m/s',color='black') #,fontsize = 6)

m.drawcoastlines(linewidth = 0.2)          
m.fillcontinents(color='gray')
plt.annotate('B',xy=(0.15,0.6),xycoords='axes fraction',fontsize=14)

print 'subplot 3'
plt.subplot(223)


z_level = 30.
height = np.zeros([JM,IM])
for i in range(0,IM):
    for j in range(0,JM):
        if h[j,i] >= -99999.:
            z_base = min(z_level,h[j,i])
            for k in range(1,KB):
                depth_k = -1. * ((h[j,i] + emean[j,i])*Z[k] + emean[j,i])
                depth_km1 = -1. * ((h[j,i] + emean[j,i])*Z[k-1] + emean[j,i])
                if depth_k <= z_base:
                    height[j,i] = height[j,i] + rhomean[k-1,j,i] * DZ[k-1] * (h[j,i] + emean[j,i]) * grav
                elif (depth_k > z_base) & (depth_km1 <= z_base):
                    ratio = (z_base - depth_km1) / (depth_k - depth_km1)
                    height[j,i] = height[j,i] + rhomean[k-1,j,i] * DZ[k-1] * (h[j,i] + emean[j,i]) * ratio * grav
            height[j,i] = height[j,i] /10000.
data = height
vmin = 0.01
vmax = 5.0
data[data<=vmin]=vmin
data[data>=vmax]=vmax
data = np.ma.masked_array(data,h<z_level)
levels = np.arange(vmin,vmax+(vmax-vmin)/30.,(vmax-vmin)/30.)
m = Basemap(llcrnrlon=lon_1,llcrnrlat=lat_1,urcrnrlon=lon_2,urcrnrlat=lat_2,\
                    resolution='h',projection='merc')
m.drawparallels(np.arange(lat_1,lat_2,stick_interval_x),labels=[1,0,0,0],linewidth=0.0)
m.drawmeridians(np.arange(lon_1,lon_2,stick_interval_y),labels=[0,0,0,1],linewidth=0.0)
plt.tick_params(axis='both',labelsize=8)
m.contourf(x,y,data,levels,cmap = cmap)
#m.colorbar()
m.contour(x,y,h,deplevels,linewidths = 0.2)
plt_zlevel = z_level
um = slice(x,y,h,ZZ,umean,plt_zlevel)
vm = slice(x,y,h,ZZ,vmean,plt_zlevel)
u_plt = np.ma.masked_array(um,mask = (h*(-1.*ZZ[-1])<plt_zlevel))
v_plt = np.ma.masked_array(vm,mask = (h*(-1.*ZZ[-1])<plt_zlevel))
xx = x[~np.isnan(u_plt)]
yy = y[~np.isnan(u_plt)]
uu = u_plt[~np.isnan(u_plt)]
vv = v_plt[~np.isnan(u_plt)]
u_plt = griddata((xx,yy),uu,(lon,lat),method='cubic')
v_plt = griddata((xx,yy),vv,(lon,lat),method='cubic')
uabs = np.sqrt(u_plt**2 + v_plt**2)
u_plt[uabs<0.001] = np.nan
v_plt[uabs<0.001] = np.nan
uv = np.sqrt(u_plt*u_plt + v_plt*v_plt)
Q = m.quiver(lon[uv>=ubig],lat[uv>=ubig],u_plt[uv>=ubig],v_plt[uv>=ubig],scale=8,color='purple')
qk = plt.quiverkey(Q,0.2,0.5,0.5,'0.5m/s',color='purple') #,fontsize = 6)
Q = m.quiver(lon[uv<ubig],lat[uv<ubig],u_plt[uv<ubig],v_plt[uv<ubig],scale=3,color='black')
qk = plt.quiverkey(Q,0.2,0.4,0.2,'0.2m/s',color='black') #,fontsize = 6)

m.drawcoastlines(linewidth = 0.2)          
m.fillcontinents(color='gray')
plt.annotate('C',xy=(0.15,0.6),xycoords='axes fraction',fontsize=14)

print 'subplot 4'
plt.subplot(224)

z_level = 50.
height = np.zeros([JM,IM])
for i in range(0,IM):
    for j in range(0,JM):
        if h[j,i] >= -99999.:
            z_base = min(z_level,h[j,i])
            for k in range(1,KB):
                depth_k = -1. * ((h[j,i] + emean[j,i])*Z[k] + emean[j,i])
                depth_km1 = -1. * ((h[j,i] + emean[j,i])*Z[k-1] + emean[j,i])
                if depth_k <= z_base:
                    height[j,i] = height[j,i] + rhomean[k-1,j,i] * DZ[k-1] * (h[j,i] + emean[j,i]) * grav
                elif (depth_k > z_base) & (depth_km1 <= z_base):
                    ratio = (z_base - depth_km1) / (depth_k - depth_km1)
                    height[j,i] = height[j,i] + rhomean[k-1,j,i] * DZ[k-1] * (h[j,i] + emean[j,i]) * ratio * grav
            height[j,i] = height[j,i] /10000.
data = height
vmin = 0.01
vmax = 5.0
data[data<=vmin]=vmin
data[data>=vmax]=vmax
data = np.ma.masked_array(data,h<z_level)
levels = np.arange(vmin,vmax+(vmax-vmin)/30.,(vmax-vmin)/30.)
m = Basemap(llcrnrlon=lon_1,llcrnrlat=lat_1,urcrnrlon=lon_2,urcrnrlat=lat_2,\
                    resolution='h',projection='merc')
m.drawparallels(np.arange(lat_1,lat_2,stick_interval_x),labels=[1,0,0,0],linewidth=0.0)
m.drawmeridians(np.arange(lon_1,lon_2,stick_interval_y),labels=[0,0,0,1],linewidth=0.0)
plt.tick_params(axis='both',labelsize=8)
m.contourf(x,y,data,levels,cmap = cmap)
# m.colorbar()
m.contour(x,y,h,deplevels,linewidths = 0.2)
plt_zlevel = z_level
um = slice(x,y,h,ZZ,umean,plt_zlevel)
vm = slice(x,y,h,ZZ,vmean,plt_zlevel)
u_plt = np.ma.masked_array(um,mask = (h*(-1.*ZZ[-1])<plt_zlevel))
v_plt = np.ma.masked_array(vm,mask = (h*(-1.*ZZ[-1])<plt_zlevel))
xx = x[~np.isnan(u_plt)]
yy = y[~np.isnan(u_plt)]
uu = u_plt[~np.isnan(u_plt)]
vv = v_plt[~np.isnan(u_plt)]
u_plt = griddata((xx,yy),uu,(lon,lat),method='cubic')
v_plt = griddata((xx,yy),vv,(lon,lat),method='cubic')
uabs = np.sqrt(u_plt**2 + v_plt**2)
u_plt[uabs<0.001] = np.nan
v_plt[uabs<0.001] = np.nan
uv = np.sqrt(u_plt*u_plt + v_plt*v_plt)
Q = m.quiver(lon[uv>=ubig],lat[uv>=ubig],u_plt[uv>=ubig],v_plt[uv>=ubig],scale=8,color='purple')
qk = plt.quiverkey(Q,0.2,0.5,0.5,'0.5m/s',color='purple') #,fontsize = 6)
Q = m.quiver(lon[uv<ubig],lat[uv<ubig],u_plt[uv<ubig],v_plt[uv<ubig],scale=3,color='black')
qk = plt.quiverkey(Q,0.2,0.4,0.2,'0.2m/s',color='black') #,fontsize = 6)

m.drawcoastlines(linewidth = 0.2)          
m.fillcontinents(color='gray')
plt.annotate('D',xy=(0.15,0.6),xycoords='axes fraction',fontsize=14)
'''
# plt.show()
figname = case + '.png'
plt.savefig(figname, dpi = 300, bbox_inches = 'tight')

plt.close()
print 'work done!'




