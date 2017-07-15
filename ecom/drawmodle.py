# -*- coding: utf-8 -*-
import numpy as np
import struct
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os
import datetime as dt
from scipy.interpolate import interp1d
from scipy.interpolate import griddata
import ResiSedFluxFieldDrawing as rsffd
import globalvar as gl


inptdict={}.fromkeys(('IN_DIRE','OUT_DIRE','IYEAR','IMONTH','IDAY0','IM','JM','KB','LOG_TSR_EL','LOG_OBS_EL','LOG_TSR_VEL','LOG_OBS_VEL','TSR_LAYER_VEL','LOG_TSR_S','LOG_OBS_S','TSR_LAYER_S','LOG_TSR_SEC','TSR_BEG',
'TSR_END','TSR_LAG','TSR_SEC_BEG','TSR_SEC_END','TSR_SEC_LAG','LOG_FIELD','N_FPT','DMIN','LOG_OPT','LOG_EPT','EPT_TIME','LOG_VPT_UV','LOG_VPT_SD','VPT_TIME','VPT_LAYER','VPT_INTERVAL','VPT_SCALE',
'LOG_SPT','SPT_TIME','SPT_LAYER','LOG_RSPT','RSPT_LAYER','LOG_SED','SED_TIME','SED_LAYER','LOG_FIELD_RCUR','RCUR_LAYER','LOG_FIELD_RSED','RSED_LAYER','RS_MIN_DEP','LOG_TAU','TAU_TIME','FPT_COORDINATE',
'FPT_XMIN','FPT_XMAX','FPT_YMIN','FPT_YMAX','FPT_COLOR_TYPE','LOG_SEC_S','SEC_S_TIME','LOG_SEC_RS','LOG_SEC_SED','SEC_SED_TIME','LOG_SEC_RSED_FLUX','SEC_NUM','SEC_CONTROL_POINTS','SEC_RESOLUTION','SEC_MAX_HIGHT'))


lndict=dict((['IN_DIRE',3],['OUT_DIRE',4],['IYEAR',6],['IMONTH',7],['IDAY0',8],['IM',10],['JM',11],['KB',12],
['LOG_TSR_EL',14],['LOG_OBS_EL',15],['LOG_TSR_VEL',16],['LOG_OBS_VEL',17],['TSR_LAYER_VEL',18],['LOG_TSR_S',19],['LOG_OBS_S',20],
['TSR_LAYER_S',21],['LOG_TSR_SEC',22],['TSR_BEG',23],['TSR_END',24],['TSR_LAG',25],['TSR_SEC_BEG',26],['TSR_SEC_END',27],['TSR_SEC_LAG',28],
['LOG_FIELD',30],['N_FPT',31],['DMIN',32],['LOG_OPT',33],['LOG_EPT',34],['EPT_TIME',35],['LOG_VPT_UV',36],['LOG_VPT_SD',37],
['VPT_TIME',38],['VPT_LAYER',39],['VPT_INTERVAL',40],['VPT_SCALE',41],['LOG_SPT',42],['SPT_TIME',43],['SPT_LAYER',44],['LOG_RSPT',45],
['RSPT_LAYER',46],['LOG_SED',47],['SED_TIME',48],['SED_LAYER',49],['LOG_FIELD_RCUR',50],['RCUR_LAYER',51],['LOG_FIELD_RSED',52],
['RSED_LAYER',53],['RS_MIN_DEP',54],['LOG_TAU',55],['TAU_TIME',56],['FPT_COORDINATE',57],['FPT_XMIN',58],['FPT_XMAX',59],['FPT_YMIN',60],
['FPT_YMAX',61],['FPT_COLOR_TYPE',62],['LOG_SEC_S',64],['SEC_S_TIME',65],['LOG_SEC_RS',66],['LOG_SEC_SED',67],['SEC_SED_TIME',68],
['LOG_SEC_RSED_FLUX',69],['SEC_NUM',70],['SEC_CONTROL_POINTS',71],['SEC_RESOLUTION',73],['SEC_MAX_HIGHT',74]))

inptnm=['IN_DIRE','OUT_DIRE','IYEAR','IMONTH','IDAY0','IM','JM','KB','LOG_TSR_EL','LOG_OBS_EL','LOG_TSR_VEL','LOG_OBS_VEL','TSR_LAYER_VEL','LOG_TSR_S','LOG_OBS_S','TSR_LAYER_S','LOG_TSR_SEC','TSR_BEG',
'TSR_END','TSR_LAG','TSR_SEC_BEG','TSR_SEC_END','TSR_SEC_LAG','LOG_FIELD','N_FPT','DMIN','LOG_OPT','LOG_EPT','EPT_TIME','LOG_VPT_UV','LOG_VPT_SD','VPT_TIME','VPT_LAYER','VPT_INTERVAL','VPT_SCALE',
'LOG_SPT','SPT_TIME','SPT_LAYER','LOG_RSPT','RSPT_LAYER','LOG_SED','SED_TIME','SED_LAYER','LOG_FIELD_RCUR','RCUR_LAYER','LOG_FIELD_RSED','RSED_LAYER','RS_MIN_DEP','LOG_TAU','TAU_TIME','FPT_COORDINATE',
'FPT_XMIN','FPT_XMAX','FPT_YMIN','FPT_YMAX','FPT_COLOR_TYPE','LOG_SEC_S','SEC_S_TIME','LOG_SEC_RS','LOG_SEC_SED','SEC_SED_TIME','LOG_SEC_RSED_FLUX','SEC_NUM','SEC_CONTROL_POINTS','SEC_RESOLUTION','SEC_MAX_HIGHT']


#============ Reading Draw Setting File =======================
with open('Draw_setting.txt') as f:
        lines = f.readlines()


for key in lndict.keys():
    data = lines[lndict[key]-1].split()
#    print 'key= %s, value= %s' % (key,data[1])
    inptdict[key]=data[1]
    
SEC_CONTROL_POINT = lines[71].split(';')    
#================= End Of Draw Setting File Reading =======================
gl.IN_DIRE=inptdict['IN_DIRE']
gl.OUT_DIRE=inptdict['OUT_DIRE'] 
gl.IYEAR=inptdict['IYEAR']		
gl.IMONTH=inptdict['IMONTH']	
gl.IDAY0=inptdict['IDAY0']		
gl.IM=inptdict['IM']	
gl.JM=inptdict['JM']	
gl.KB=inptdict['KB']
gl.LOG_TSR_EL=inptdict['LOG_TSR_EL']	
gl.LOG_OBS_EL=inptdict['LOG_OBS_EL']
gl.LOG_TSR_VEL=inptdict['LOG_TSR_VEL']
gl.LOG_OBS_VEL=inptdict['LOG_OBS_VEL']
gl.TSR_LAYER_VEL=inptdict['TSR_LAYER_VEL']
gl.LOG_TSR_S=inptdict['LOG_TSR_S']
gl.LOG_OBS_S=inptdict['LOG_OBS_S']
gl.TSR_LAYER_S=inptdict['TSR_LAYER_S'] 
gl.LOG_TSR_SEC=inptdict['LOG_TSR_SEC']
gl.TSR_BEG=inptdict['TSR_BEG']
gl.TSR_END=inptdict['TSR_END']
gl.TSR_LAG=inptdict['TSR_LAG']
gl.TSR_SEC_BEG=inptdict['TSR_SEC_BEG']
gl.TSR_SEC_END=inptdict['TSR_SEC_END']
gl.TSR_SEC_LAG=inptdict['TSR_SEC_LAG']
gl.LOG_FIELD=inptdict['LOG_FIELD']
gl.N_FPT=inptdict['N_FPT']
gl.DMIN=inptdict['DMIN']
gl.LOG_OPT=inptdict['LOG_OPT']
gl.LOG_EPT=inptdict['LOG_EPT']
gl.EPT_TIME=inptdict['EPT_TIME']
gl.LOG_VPT_UV=inptdict['LOG_VPT_UV']
gl.LOG_VPT_SD=inptdict['LOG_VPT_SD']
gl.VPT_TIME=inptdict['VPT_TIME']
gl.VPT_LAYER=inptdict['VPT_LAYER']
gl.VPT_INTERVAL=inptdict['VPT_INTERVAL']
gl.VPT_SCALE=inptdict['VPT_SCALE']
gl.LOG_SPT=inptdict['LOG_SPT']
gl.SPT_TIME=inptdict['SPT_TIME']
gl.SPT_LAYER=inptdict['SPT_LAYER']
gl.LOG_RSPT=inptdict['LOG_RSPT']
gl.RSPT_LAYER=inptdict['RSPT_LAYER']
gl.LOG_SED=inptdict['LOG_SED']
gl.SED_TIME=inptdict['SED_TIME']
gl.SED_LAYER=inptdict['SED_LAYER']
gl.LOG_FIELD_RCUR=inptdict['LOG_FIELD_RCUR']
gl.RCUR_LAYER=inptdict['RCUR_LAYER']
gl.LOG_FIELD_RSED=inptdict['LOG_FIELD_RSED']
gl.RSED_LAYER=inptdict['RSED_LAYER']
gl.RS_MIN_DEP=inptdict['RS_MIN_DEP']
gl.LOG_TAU=inptdict['LOG_TAU']
gl.TAU_TIME=inptdict['TAU_TIME']
gl.FPT_COORDINATE=inptdict['FPT_COORDINATE']
gl.FPT_XMIN=inptdict['FPT_XMIN']
gl.FPT_XMAX=inptdict['FPT_XMAX']
gl.FPT_YMIN=inptdict['FPT_YMIN']
gl.FPT_YMAX=inptdict['FPT_YMAX']
gl.FPT_COLOR_TYPE=inptdict['FPT_COLOR_TYPE']
gl.LOG_SEC_S=inptdict['LOG_SEC_S']
gl.SEC_S_TIME=inptdict['SEC_S_TIME']
gl.LOG_SEC_RS=inptdict['LOG_SEC_RS']
gl.LOG_SEC_SED=inptdict['LOG_SEC_SED']
gl.SEC_SED_TIME=inptdict['SEC_SED_TIME']
gl.LOG_SEC_RSED_FLUX=inptdict['LOG_SEC_RSED_FLUX']
gl.SEC_NUM=inptdict['SEC_NUM']
gl.SEC_CONTROL_POINTS=inptdict['SEC_CONTROL_POINTS']
gl.SEC_RESOLUTION=inptdict['SEC_RESOLUTION']
gl.SEC_MAX_HIGHT=inptdict['SEC_MAX_HIGHT']

#================== Time Seriers Drawing ====================
#IM=int(gl.IM)
#JM=int(gl.JM)
#KB=int(gl.KB)

case = 'case2011'
indir = './' + case + '/Ecomsi3dsed/output/field_distri/'
kind1 = "temperature/t_field_"
kind2 = 'salinity/s_field_'
kind3 = 'current/v_field_'
kind4 = 'elevation/el_field_'
kind5 = 'momentum/mmt_field_'

ntime1 = 8761
ntime2 = ntime1+24
ntime = ntime2-ntime1+1
factor = 1.0e5
date0 = dt.datetime(2008,1,1)
grav = 9.801

date = date0 + dt.timedelta(hours=ntime1)
tracer = 0
total_tracers = 2
pressure = 0.
lon_1 = 119
lon_2 = 126
lat_1 = 25
lat_2 = 33
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


#  read grid file
with open('./grid_info/ch_hzbc.grd') as f:
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
with open('./grid_info/ch_hzbc_dep.dat') as f:
    deplines = f.readlines()

startline = 1
for ii in range(startline,startline+IMJM):
    data = deplines[ii].split()
    i = int(data[0])-1
    j = int(data[1])-1
    h[j,i] = data[2]    


#    hr=f.read(IM*JM); #读水深数据
#    xr=f.read(IM*JM); #读网格中心点的X坐标
#    yr=f.read(IM*JM); #读网格中心点的Y坐标
#    erecord=f.read(1); #读记录大小信息
#    f.close();


    
#if (gl.FPT_COLOR_TYPE=='C'): #按照绘图的颜色类型设置底图的背景色
#   bgcolor=[1,1,1]
#elif(gl.FPT_COLOR_TYPE=='G'):
#   bgcolor=[0.3,0.3,0.3]

     
#hr =np.reshape(hr,IM,JM); #转换成相应的矩阵形式
#xr=np.reshape(xr,IM,JM); #转换成相应的矩阵形式
#yr=np.reshape(yr,IM,JM); #转换成相应的矩阵形式
    

#xr=xr/1000; #将54坐标单位从m转换为km
#yr=yr/1000; #将54坐标单位从m转换为km
'''
#================== reade grid data ====================




def TimeSeriersDrawing():
    import globalvar as gl    
    if (gl.LOG_TSR_EL=='T' & gl.LOG_TSR_VEL=='T'& gl.LOG_TSR_S=='T'& gl.LOG_TSR_SEC=='T'):
       pass
#================== End Of Time Seriers Drawing ====================
    
#================== Section Time Seriers Drawing =================== 
#================== End Of Section Time Seriers Drawing ================ 


#=============== Field and Section Distribution Drawing =================

def  FieldDistributionDrawing():
    import globalvar as gl 
    if ( gl.LOG_FIELD =='T'): #平面断面绘图开关（T画，F不画）
        pass
   
#------------- Salinity Section Distribution Drawing -------------------

def  SalSecDistribtionDrawing():
    import globalvar as gl 
    if (gl.LOG_SEC_S=='T'): #如果需要画盐度的断面图
       pass
 #------------- End Salinity Section Distribution Drawing ----------------
 
 
 #------------- Residual Salinity Section Distribution Drawing -----------
   #如果需要画盐度的断面图
def ResiSalSecDrawing():
    import globalvar as gl 
    if (gl.LOG_SEC_RS=='T'):
        pass 
 #------------- End Residual Salinity Section Distribution Drawing -------- 
    
    
#------------- Sediment Section Distribution Drawing -------------------

def SedSecDistribtionDrawing():
    import globalvar as gl 
    if (gl.LOG_SEC_SED =='T'): #如果需要画含沙量的断面图      
        pass
#------------- End Sediment Section Distribution Drawing ---------------- 
       
#------------- Residual Sediment Flux Field Distribution Drawing -----------


      

def ResiCurFluxFieldDrawing():  
  if (inptdict['LOG_FIELD_RSED']=='T'): #如果需要画余含沙量断面通量图
     print '===== Residual Sediment Flux Field Distribution Drawing ======'

  
#------------- End Residual Sediment Flux Field Distribution Drawing ----------------  
 
#------------- Residual Sediment Flux Field Distribution Drawing -----------
def ResiSedFluxFieldDrawing():  
#   import globalvar as gl
    import ResiSedFluxFieldDrawing as ResiSedF
    print '===== Run ResiSedFluxFieldDrawing in main1 ======'
    plt1=ResiSedF.ResiSedFluxFieldDrawing()
    plt1.draw_ResiSedFluxField()



 
      nisland=13;
    eval(['load LandFiles_XY\land.dat']);
    for k1 = 1:nisland;
        eval(['load LandFiles_XY\island' int2str(k1) '.dat']);
    end
    
    eval(['load LandFiles_XY\prj\nanhuibiantan_weiken.dat']);
    
    eval(['load LandFiles_XY\prj\dnc_2phase.dat']);
    eval(['load LandFiles_XY\prj\yh101.dat']);

        land=land/1000;
        
        for k1 = 1:nisland;
            eval(['island' int2str(k1) '=island' int2str(k1) '/1000;']);
        end
        
        nanhuibiantan_weiken=nanhuibiantan_weiken/1000;
        dnc_2phase=dnc_2phase/1000;
        yh101=yh101/1000;
#----------------------END------------------------
    
    print '===== Run ResiSedFluxFieldDrawing in main2 ======'
    
#    fid_XYH_BIN=fopen('ch_hzbc_griddepth'); #打开网格中心坐标、水深文件，该文件由模式输出

#------------- End Residual Sediment Flux Field Distribution Drawing ----------------  

#------------- Residual Sediment Flux Section Distribution Drawing -----------


def ResiSedFluxSecDrawing(): 
    import globalvar as gl 
    print '===== End  Residual Sediment Flux Section Distribution Drawing ======'
    if (gl.LOG_SEC_RSED_FLUX=='T'): #如果需要画余含沙量断面通量图
       pass

#------------- End Residual Sediment Flux Section Distribution Drawing ----------------    

if __name__ == '__main__':
    if (gl.LOG_FIELD_RSED=='T'): #如果需要画余含沙量断面通量图
        print '=====  Residual Sediment Flux Field Distribution Drawing in __main__======'
        pass
#        ResiSedFluxFieldDrawing()
'''    