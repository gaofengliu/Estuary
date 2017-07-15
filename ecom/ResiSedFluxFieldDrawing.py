# -*- coding: utf-8 -*-
import globalvar as gl


class ResiSedFluxFieldDrawing(object):
    def __init__(self):
            self.IYEAR=0
            self.IMONTH=0
            self.IDAY0=0
#  if (LOG_FIELD_RSED=='T') #如果需要画余含沙量断面通量图
    #----------读取网格中心坐标和水深数据-------------------------
            self.Init_day=0 #计算模式起算的天数 
            
#Init_day=datenum([IYEAR IMONTH IDAY0])
    def draw_ResiSedFluxField(self):
        print 'Run draw_ResiSedFluxField in class'
        pass
''' 
#   fid_XYH_BIN=fopen('ch_hzbc_griddepth'); #打开网格中心坐标、水深文件，该文件由模式输出
    brecord=fread(fid_XYH_BIN,1,'integer*4'); #读记录大小信息
    h=fread(fid_XYH_BIN,IM*JM,'real*4'); #读水深数据
    xr=fread(fid_XYH_BIN,IM*JM,'real*4'); #读网格中心点的X坐标
    yr=fread(fid_XYH_BIN,IM*JM,'real*4'); #读网格中心点的Y坐标
    erecord=fread(fid_XYH_BIN,1,'integer*4'); #读记录大小信息
    fclose(fid_XYH_BIN);
    
    if (brecord~=erecord) #如果一条记录的前后所记录的大小不一致，则退出程序
        disp('Error in reading file ch_hzbc_griddepth, please check!')
#       break

    
   if (FPT_COLOR_TYPE=='C') #按照绘图的颜色类型设置底图的背景色
        bgcolor=[1 1 1];
    elseif (FPT_COLOR_TYPE=='G')
        bgcolor=[0.3 0.3 0.3];

     
    h=reshape(h,IM,JM); #转换成相应的矩阵形式
    xr=reshape(xr,IM,JM); #转换成相应的矩阵形式
    yr=reshape(yr,IM,JM); #转换成相应的矩阵形式
    

        xr=xr/1000; #将54坐标单位从m转换为km
        yr=yr/1000; #将54坐标单位从m转换为km
    
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
        disp('Residual Sediment Flux Field Distribution Drawing ...')
        RSED_Flux_Sec_Pic_Path=fullfile(OUT_DIRE,'residual_distri\field\sediment_flux'); #余含沙量断面通量图存放路径（指定）
        [s,mess,messid]=mkdir(RSED_Flux_Sec_Pic_Path); #建立余含沙量断面图的存放文件夹（指定）
        
        REPT_Comp_Filepath=fullfile(IN_DIRE,'resi_euler'); #模式的余水位场计算结果文件路径
        REPT_Comp_Fileinfo=dir(fullfile(REPT_Comp_Filepath,'resi_euler_current_*.out')); #模式计算的余水位场（包括文件名，修改时间，大小，是否为文件路径）
        RSED_Flux_Comp_Filepath=fullfile(IN_DIRE,'resi_flux'); #模式的余含沙量断面通量计算结果文件路径
        RSED_Flux_Comp_Fileinfo=dir(fullfile(RSED_Flux_Comp_Filepath,'resi_sed_flux_3d_*.out')); #模式计算的含沙量断面通量（包括文件名，修改时间，大小，是否为文件路径）        
        RSED_Flux_contourscale=[-0.1 0 0.1 0.5 1 1.5 2 5]; #需要画的等值线        
        #计算网格边长
        
        for i=2:IM
            for j=2:JM
                h1(i,j)=sqrt((xr(i-1,j)-xr(i,j))^2+(yr(i-1,j)-yr(i,j))^2);
                h2(i,j)=sqrt((xr(i,j-1)-xr(i,j))^2+(yr(i,j-1)-yr(i,j))^2);
            end
        end
       
        for i=1:length(RSED_Flux_Comp_Fileinfo) #按照需要画图的数量进行循环
            
            #-----------------Residual Elevation Field Distribution Data Reading ------------------
            REPT_Comp_Filename=REPT_Comp_Fileinfo(i).name; #模式计算输出的余水位文件  
            eval(['fid_REPT_Comp = fopen(' fullfile(REPT_Comp_Filepath,REPT_Comp_Filename) ',''r'');']); #打开余水位文件
            REPT_Comp_Data=nan*ones(IM,JM);
            Note=fscanf(fid_REPT_Comp,'#s',7);
            REPT_B=fscanf(fid_REPT_Comp,'%f',1);
            Note=fscanf(fid_REPT_Comp,'#s',2);
            REPT_E=fscanf(fid_REPT_Comp,'%f',1);
            Note=fgetl(fid_REPT_Comp);
            Note=fgetl(fid_REPT_Comp);         
            
            while(~feof(fid_REPT_Comp)) #判读是否读到文件尾部
                I=fscanf(fid_REPT_Comp,'%d',1); #读I
                J=fscanf(fid_REPT_Comp,'%d',1); #读J
                Note=fscanf(fid_REPT_Comp,'%f',2); #读X,Y
                REPT_Comp_Data(I,J)=fscanf(fid_REPT_Comp,'%f',1); #读X,Y位置的水位值
                Note=fgetl(fid_REPT_Comp); #读掉其余部分
            end
            
            fclose(fid_REPT_Comp); #关闭计算输出的余水位文件
            #------------------- End Residal Elevation Field Distribution Data Reading -----------------------
           
           
            #-----------------Residual Sediment flux Field Distribution Data Reading ------------------
            RSED_Flux_Comp_Filename=RSED_Flux_Comp_Fileinfo(i).name; #模式计算输出的余含沙量文件  
            eval(['fid_RSED_Flux_Comp = fopen(' fullfile(RSED_Flux_Comp_Filepath,RSED_Flux_Comp_Filename) ',''r'');']); #打开余含沙量文件
            RSED_Flux_U_Comp_Data=nan*ones(IM,JM,KB-1);
            RSED_Flux_V_Comp_Data=nan*ones(IM,JM,KB-1);
            Note=fscanf(fid_RSED_Flux_Comp,'#s',7);
            RSED_Flux_B=fscanf(fid_RSED_Flux_Comp,'%f',1);
            Note=fscanf(fid_RSED_Flux_Comp,'#s',2);
            RSED_Flux_E=fscanf(fid_RSED_Flux_Comp,'%f',1);
            Note=fgetl(fid_RSED_Flux_Comp);
            Note=fgetl(fid_RSED_Flux_Comp);
            disp(['  Residual Sediment Section Distribution from ' num2str(RSED_Flux_B) ' to ' num2str(RSED_Flux_E) ' Hour Drawing'])
            
            while(~feof(fid_RSED_Flux_Comp)) #判读是否读到文件尾部
                I=fscanf(fid_RSED_Flux_Comp,'%d',1); #读I
                J=fscanf(fid_RSED_Flux_Comp,'%d',1); #读J
                Note=fscanf(fid_RSED_Flux_Comp,'%f',2+KB-1); #读X,Y和余含沙量数据
                for kk=1:KB-1
                    RSED_Flux_U_Comp_Data(I,J,kk)=fscanf(fid_RSED_Flux_Comp,'%f',1); #读X,Y位置U方向的余含沙通量
                    RSED_Flux_V_Comp_Data(I,J,1:kk)=fscanf(fid_RSED_Flux_Comp,'%f',1); #读X,Y位置V方向的余含沙通量
                end
                Note=fgetl(fid_RSED_Flux_Comp); #读掉末尾的回车
            end
            
            fclose(fid_RSED_Flux_Comp); #关闭计算输出的余含沙量文件

            #------------------- End Residual Sediment flux Field Distribution Data Reading -----------------------
            
            fsm=ones(IM,JM);
            for i1=1:IM
                for j1=1:JM
                    if isnan(RSED_Flux_U_Comp_Data(i1,j1))==1 #按照fsmadd判断
                        fsm(i1,j1)=0; 
                        REPT_Comp_Data(i1,j1)=-SEC_MAX_HIGHT; 
                        RSED_Flux_U_Comp_Data(i1,j1,:)=NaN;    
                        RSED_Flux_V_Comp_Data(i1,j1,:)=NaN;    
                    end
                    
                    if (h(i1,j1)<RS_MIN_DEP) #小于RS_MIN_DEP位置不画
                        fsm(i1,j1)=0;
                        REPT_Comp_Data(i1,j1)=-SEC_MAX_HIGHT;
                        RSED_Flux_U_Comp_Data(i1,j1,:)=NaN;    
                        RSED_Flux_V_Comp_Data(i1,j1,:)=NaN;    
                    end
                end
            end  
        
            Pic_time_b=datevec(Init_day+RSED_Flux_B/24); #根据模式设置的起始时间计算余含沙量场统计的起始时间
            Pic_time_e=datevec(Init_day+RSED_Flux_E/24); #根据模式设置的起始时间计算余含沙量场统计的结束时间
            

            
      #=================== 绘制余沙场分布图 ======================
      
        NRSED_LAYER=str2num(RSED_LAYER); #将读入的含沙量绘图层数的字符数组转化为浮点数组
        Num_RSED_LAYER=length(NRSED_LAYER); #总共需要画几个层次的图像
                   
            for k=1:Num_RSED_LAYER
                close(figure(1));
                disp(['    Layer ' num2str(NRSED_LAYER(k))])
                
            figure(1),clf;
            if (FPT_COLOR_TYPE=='C') #根据绘图的颜色类型，选择颜色包文件
                map=load ('Colormap\cm_sed_C.dat');
            elseif (FPT_COLOR_TYPE=='G')
                map=load ('Colormap\cm_sed_G.dat');
            end
            map=[bgcolor;map]; #添加背景颜色，让含沙量为nan的位置按背景色画
            colormap(map);  #读入颜色包
           
#             scrz = get(0,'ScreenSize');
#             SPPI = get(0,'ScreenPixelsPerInch');
#             #绘制长12cm，高8cm图像的像素
#             W = SPPI*22/2.54;
#             H = SPPI*24/2.54;
#             X0 = (scrz(3)-W)/2;#绘制在屏幕中间
#             Y0 = (scrz(4)-H)/2;
#             set(gcf,'Position',[X0 Y0 W H]);
            
           set(gcf,'position',[50 50 1100 600]);
            x0=0.1;  y0=0.1; width=0.85; height=0.8;                    
           subplot('position',[x0 y0 width height]);
            axis equal;
            
            set(gca,'box','on','FontName','times new roman','FontSize',12);
            set(gca,'color',bgcolor);
          
            set(gcf,'inverthardcopy','off'); #保存图片时按照设置的颜色而不自动调节
            set(gcf,'color',[1 1 1]);        
            hold on; 
            
            
             set(gca,'xlim',[FPT_XMIN FPT_XMAX],'xtick',linspace(FPT_XMIN,FPT_XMAX,4));
             set(gca,'ylim',[FPT_YMIN FPT_YMAX],'ytick',linspace(FPT_YMIN,FPT_YMAX,4));   
             
                                     
            warning off;
#    画contourf时会有很多警告信息，但不影响画图，所以关掉warning
#           contourf(xr,yr,sqrt(RSED_Flux_U_Comp_Data(:,:,NRSED_LAYER(k)).^2+RSED_Flux_V_Comp_Data(:,:,NRSED_LAYER(k)).^2),RSED_Flux_contourscale);
#             shading flat; 
#             caxis([-0.1 5]); #为了把背景色放到颜色包，而又不显示在colorbar里，‘/20’是由于颜色包里有20个颜色
#             if (FPT_COLOR_TYPE=='C') #如果画彩色图，则需要画colorbar
#                 hc=colorbar;
#                 set(hc,'ylim',[0 5],'FontName','times new roman','FontSize',12); #设置colorbar的显示范围
#             end            
#            [cc,hh]=contour(xr,yr,sqrt(RSED_Flux_U_Comp_Data(:,:,NRSED_LAYER(k)).^2+RSED_Flux_V_Comp_Data(:,:,NRSED_LAYER(k)).^2),RSED_Flux_contourscale);
#            clabel(cc,hh,'FontSize',10,'FontName','Times new roman');
#            set(hh,'EdgeColor',[0.3 0.3 0.3],'linewidth',1) #设置等值线的颜色和粗细
            
                if (NRSED_LAYER(k)==1)
                    RSED_LAYER_NAME='Surface';
                elseif (NRSED_LAYER(k)==KB-1)
                    RSED_LAYER_NAME='Bottom';
                elseif (NRSED_LAYER(k)==2)
                    RSED_LAYER_NAME='2nd layer';
                elseif (NRSED_LAYER(k)==3)
                    RSED_LAYER_NAME='3rd layer';
                else
                    RSED_LAYER_NAME=[num2str(NRSED_LAYER(k)) 'th layer'];
                end
            
      #     title(['Residual sediment flux during time from ' num2str(Pic_time_b(4)) ':00 ' num2str(Pic_time_b(2),:) ' ' num2str(Pic_time_b(3)) ' ' num2str(Pic_time_b(1)) ' to ' num2str(Pic_time_e(4)) ':00 ' num2str(Pic_time_e(2),:) ' ' num2str(Pic_time_e(3)) ' ' num2str(Pic_time_e(1)) ' ' RSED_LAYER_NAME])
            hold on  

            #在quiver函数中设定箭头按实际大小画，而不是自动调节大小，所以需要把将原始数据放大到和坐标轴相应的大小尺度
            Draw_Flux_U=RSED_Flux_U_Comp_Data*(FPT_XMAX-FPT_XMIN)*0.05/VPT_SCALE; #(FPT_XMAX-FPT_XMIN)*0.05/VPT_SCALE表示将x轴的5#的长度对应VPT_SCALE的大小
            Draw_Flux_V=RSED_Flux_V_Comp_Data*(FPT_YMAX-FPT_YMIN)*0.05/VPT_SCALE; #(FPT_YMAX-FPT_YMIN)*0.05/VPT_SCALE表示将y轴的5#的长度对应VPT_SCALE的大小           
#           quiver(xr(1:VPT_INTERVAL:end,1:VPT_INTERVAL:end),yr(1:VPT_INTERVAL:end,1:VPT_INTERVAL:end),Draw_Flux_U(1:VPT_INTERVAL:end,1:VPT_INTERVAL:end,1),Draw_Flux_V(1:VPT_INTERVAL:end,1:VPT_INTERVAL:end,1),0,'k');
       
            plotvector(xr(1:VPT_INTERVAL:end,1:VPT_INTERVAL:end),yr(1:VPT_INTERVAL:end,1:VPT_INTERVAL:end),Draw_Flux_U(1:VPT_INTERVAL:end,1:VPT_INTERVAL:end,NRSED_LAYER(k)),Draw_Flux_V(1:VPT_INTERVAL:end,1:VPT_INTERVAL:end,NRSED_LAYER(k)),VPT_SCALE,1,'k','line','uv','fix',1);            
            hold on
       
 
            
            #---------绘制岸线、岛------------------------
        #    fill(land(:,1),land(:,2),bgcolor); 
            plot(land(:,1),land(:,2),'color','k');

            hold on
            
            for k2 =1:nisland;
           #     eval(['fill( island' int2str(k2) '(:,1),island' int2str(k2) '(:,2),[' num2str(bgcolor) ']);']);
                eval(['plot( island' int2str(k2) '(:,1),island' int2str(k2) '(:,2),''color'',''k'');']);
                hold on
            end         
            
         #   fill(nanhuibiantan_weiken(:,1),nanhuibiantan_weiken(:,2),bgcolor);
            plot(nanhuibiantan_weiken(:,1),nanhuibiantan_weiken(:,2),'color','k');
            hold on
            
            plot(dnc_2phase(:,1),dnc_2phase(:,2),'color','k');  
            hold on
  
            plot(yh101(:,1),yh101(:,2),'color','k');  
            hold on
            
            x_scale=FPT_XMIN+(FPT_XMAX-FPT_XMIN)*0.05; #横向为从左向右5#的位置
            y_scale=FPT_YMIN+(FPT_YMAX-FPT_YMIN)*0.12; #纵向为从下而上12#的位置

            quiver([x_scale x_scale],[y_scale y_scale],[(FPT_XMAX-FPT_XMIN)*0.05 0],[0 (FPT_YMAX-FPT_YMIN)*0.05],0,'k'); #按照画箭头尺度（这里画图片纵横各5#）
          
            hold on
            text(FPT_XMIN+(FPT_XMAX-FPT_XMIN)*0.05,FPT_YMIN+(FPT_YMAX-FPT_YMIN)*0.075,[num2str(VPT_SCALE) 'kg/m^2s'],'FontName','times new roman','FontSize',10); #标注箭头的量值
            hold on
            
#              text(FPT_XMIN+(FPT_XMAX-FPT_XMIN)*0.352,FPT_YMIN+(FPT_YMAX-FPT_YMIN)*0.41,'九段沙','FontSize',10); 
#              text(FPT_XMIN+(FPT_XMAX-FPT_XMIN)*0.135,FPT_YMIN+(FPT_YMAX-FPT_YMIN)*0.92,'横沙岛','FontSize',10);
#              hold on
            

              if (FPT_COORDINATE=='BL') #按照图像的坐标轴选取绘制相应的label
              
             [FPTBL_XMIN,FPTBL_YMIN]=xy2bl(123,FPT_XMIN*1000,FPT_YMIN*1000);
             [FPTBL_XMAX,FPTBL_YMAX]=xy2bl(123,FPT_XMAX*1000,FPT_YMAX*1000); 
             
             formatSpec = '#10.2f\n';
             xlon=num2str(linspace(FPTBL_XMIN,FPTBL_XMAX,4)',formatSpec);
             ylat=num2str(linspace(FPTBL_YMIN,FPTBL_YMAX,4)',formatSpec);

             set(gca,'xTickLabel',xlon);
             set(gca,'yTickLabel',ylat);             

             xlabel('Longitude (\circE)');
             ylabel('Latitude (\circN)');
                
            elseif (FPT_COORDINATE=='XY')                
            set(gca,'xtick',fix(linspace(FPT_XMIN,FPT_XMAX,4)));
            set(gca,'ytick',fix(linspace(FPT_YMIN,FPT_YMAX,4)));
                xlabel('Distance (km)');
                ylabel('Distance (km)');
              end
             
            
             figure(1),eval(['print(gcf,''-dpng'',' fullfile(RSED_Flux_Sec_Pic_Path,RSED_Flux_Comp_Filename(1:end-4)) '_' RSED_LAYER_NAME,');']); #保存平面图像
             figure(1),eval(['saveas(gcf,'' fullfile(RSED_Flux_Sec_Pic_Path,RSED_Flux_Comp_Filename(1:end-4)) '_' RSED_LAYER_NAME '.eps'',''psc2 '');']);  #保存断面图像
#              figure(1),eval(['saveas(gcf,'' fullfile(RSED_Flux_Sec_Pic_Path,RSED_Flux_Comp_Filename(1:end-4)) '_' RSED_LAYER_NAME '.emf'');']);  #保存断面图像
                  
            end
            close(figure(1));

            #=================== 结束余沙场分布图绘制 ===============
        end   
 end
 '''