function [ output_args ] = FieldDistributionDrawing( input_args )
% %   UNTITLED3 Summary of this function goes here
% %   Detailed explanation goes here
global IN_DIRE OUT_DIRE IYEAR IMONTH IDAY0 IM JM KB;
global LOG_TSR_EL LOG_OBS_EL LOG_TSR_VEL LOG_OBS_VEL TSR_LAYER_VEL;
global LOG_TSR_S LOG_OBS_S TSR_LAYER_S LOG_TSR_SEC TSR_BEG TSR_END TSR_LAG;
global TSR_SEC_BEG TSR_SEC_END TSR_SEC_LAG N_FPT DMIN;
global LOG_OPT LOG_EPT; 
global EPT_TIME LOG_VPT_UV LOG_VPT_SD VPT_TIME VPT_LAYER VPT_INTERVAL;
global VPT_SCALE 
global LOG_FIELD LOG_SPT SPT_TIME SPT_LAYER LOG_RSPT RSPT_LAYER LOG_SED;
global SED_TIME SED_LAYER LOG_RSED RSED_LAYER RS_MIN_DEP LOG_TAU TAU_TIME;
global FPT_COORDINATE FPT_XMIN FPT_XMAX FPT_YMIN FPT_YMAX FPT_COLOR_TYPE;
global LOG_SEC_S SEC_S_TIME LOG_SEC_RS LOG_SEC_SED SEC_SED_TIME LOG_SEC_RSED_FLUX;
global SEC_P_TYPE SEC_NUM SEC_RESOLUTION SEC_MAX_HIGHT SEC_CONTROL_POINTS;
global xr yr ; 

disp('=============== Field and Section Distribution Drawing =================')

    Init_day=datenum([IYEAR IMONTH IDAY0]); %计算模式起算的天数    
    %----------读取网格中心坐标和水深数据-------------------------
    fid_XYH_BIN=fopen('ch_hzbc_griddepth'); %打开网格中心坐标、水深文件，该文件由模式输出
    brecord=fread(fid_XYH_BIN,1,'integer*4'); %读记录大小信息
    h=fread(fid_XYH_BIN,IM*JM,'real*4'); %读水深数据
    xr=fread(fid_XYH_BIN,IM*JM,'real*4'); %读网格中心点的X坐标
    yr=fread(fid_XYH_BIN,IM*JM,'real*4'); %读网格中心点的Y坐标
    erecord=fread(fid_XYH_BIN,1,'integer*4'); %读记录大小信息
    fclose(fid_XYH_BIN);
    
    if (brecord~=erecord) %如果一条记录的前后所记录的大小不一致，则退出程序
        disp('Error in reading file ch_hzbc_griddepth, please check!')
%       break
    end
    
    
    h=reshape(h,IM,JM); %转换成相应的矩阵形式
    xr=reshape(xr,IM,JM); %转换成相应的矩阵形式
    yr=reshape(yr,IM,JM); %转换成相应的矩阵形式
    
    
    fid_XYH_XY=fopen('ch_hzbc_xy_h.txt','w'); %打开新文件
    for j=1:JM
        for i=1:IM
            fprintf(fid_XYH_XY,'%20.6f %19.6f %19.6f\r\n',xr(i,j),yr(i,j),h(i,j)); %按照xr(i,j),yr(i,j),h(i,j)的格式写入网格信息
        end
    end
    fclose(fid_XYH_XY);
    
    if (FPT_COORDINATE=='XY') %如果平面图形为54坐标
        xr=xr/1000; %将54坐标单位从m转换为km
        yr=yr/1000; %将54坐标单位从m转换为km
    end
    
    if (FPT_COORDINATE=='BL') %如果平面图形需要以经纬度为坐标
        disp('Please Trans  ''ch_hzbc_xy_h.txt''  To  ''ch_hzbc_xy_h.dat''  With  ''Coor.exe''') %用坐标转化工具进行坐标转换
        pause
        
        
        fid_XYH_BL=fopen('ch_hzbc_xy_h.dat','r');
        XYH_BL=fscanf(fid_XYH_BL,'%f %f %f',[3 inf]); %读入转化好的经纬度坐标和水深数据
        xr=reshape(XYH_BL(1,:),IM,JM); %将经度写成矩阵形式
        yr=reshape(XYH_BL(2,:),IM,JM); %将纬度写成矩阵形式       
        fclose(fid_XYH_BL);
    end
    %----------------------END------------------------
    
    
    %----------读取绘图所需的岸线资料-------------------------
    nisland=13;
    eval(['load LandFiles_' FPT_COORDINATE  '\land.dat']);
    for k1 = 1:nisland;
        eval(['load LandFiles_' FPT_COORDINATE '\island' int2str(k1) '.dat']);
    end
    
    eval(['load LandFiles_' FPT_COORDINATE  '\nanhuibiantan_weiken.dat']);
    
    eval(['load LandFiles_' FPT_COORDINATE  '\deep_water_way.dat']);
    
    if (FPT_COORDINATE=='XY') %如果使用54坐标，将岸线坐标单位从m转换为km
        land=land/1000;
        
        for k1 = 1:nisland;
            eval(['island' int2str(k1) '=island' int2str(k1) '/1000;']);
        end
        
        nanhuibiantan_weiken=nanhuibiantan_weiken/1000;
        deep_water_way=deep_water_way/1000;
    end
    
    %-----------月份的英文简写---------------------
    mon(1,:)='Jan.';
    mon(2,:)='Feb.';
    mon(3,:)='Mar.';
    mon(4,:)='Apr.';
    mon(5,:)='May.';
    mon(6,:)='Jun.';
    mon(7,:)='Jul.';
    mon(8,:)='Aug.';
    mon(9,:)='Sep.';
    mon(10,:)='Oct.';
    mon(11,:)='Nov.';
    mon(12,:)='Dec.';    
    
    %----------------------END------------------------
    
    
    
    if (FPT_COLOR_TYPE=='C') %按照绘图的颜色类型设置底图的背景色
        bgcolor=[0 1 0];
    elseif (FPT_COLOR_TYPE=='G')
        bgcolor=[0.3 0.3 0.3];
    end
    
    %------------- Output Site and Section Position Drawing ---------------
    if (LOG_OPT=='T') %如果需要画输出站点和断面示意图
        disp('Output Site and Section Position Drawing ...')
        close(figure(1))
        OPT_Pic_Path=fullfile(OUT_DIRE,'output'); %输出站点和断面示意图存放路径
        [s,mess,messid]=mkdir(OPT_Pic_Path); %输出站点和断面示意图存放的文件夹
        OPT_Comp_Site_Filepath=fullfile(IN_DIRE,'timeseries'); %模式的站点计算结果文件路径
        OPT_Comp_Section_Filepath=fullfile(IN_DIRE,'secflux'); %模式的断面计算结果文件路径
        
        OPT_Comp_Site_Fileinfo=dir(fullfile(OPT_Comp_Site_Filepath,'*.out')); %模式计算站点文件信息（包括文件名，修改时间，大小，是否为文件路径）
        OPT_Comp_Section_Fileinfo=dir(fullfile(OPT_Comp_Section_Filepath,'*.out')); %模式计算断面文件信息（包括文件名，修改时间，大小，是否为文件路径）
        OPT_Site_Num=length(OPT_Comp_Site_Fileinfo); %计算有几个输出点        
        OPT_Section_Num=length(OPT_Comp_Section_Fileinfo); %计算有几个输出断面    
        
      for i=1:OPT_Section_Num
            OPT_Comp_Section_Filename=OPT_Comp_Section_Fileinfo(i).name; %计算结果文件的文件名
            eval(['fid_OPT_Comp_Section = fopen(''' fullfile(OPT_Comp_Section_Filepath,OPT_Comp_Section_Filename) ''');']); %打开站点结果文件
            
            
            note=fscanf(fid_OPT_Comp_Section,'%s',1);
            secname=fscanf(fid_OPT_Comp_Section,'%s',1);
            note=fgetl(fid_OPT_Comp_Section);
            note=fgetl(fid_OPT_Comp_Section);
            note=fscanf(fid_OPT_Comp_Section,'%s',1);
            note=fscanf(fid_OPT_Comp_Section,'%s',1);
            i1=fscanf(fid_OPT_Comp_Section,'%f',1);
            j1=fscanf(fid_OPT_Comp_Section,'%f',1);
            note=fscanf(fid_OPT_Comp_Section,'%s',1);
            i2=fscanf(fid_OPT_Comp_Section,'%f',1);
            j2=fscanf(fid_OPT_Comp_Section,'%f',1);
            
            fclose(fid_OPT_Comp_Section); %关闭模式计算结果文件
            
            plot([xr(i1,j1) xr(i2,j2)],[yr(i1,j1),yr(i2,j2)],'linewidth',1.5,'color','k');
            
            hold on
            text(xr(i1,j1),yr(i1,j1),secname(5:end),'FontName','times new roman','FontSize',12)
            
        end
        
        
        
        
        set(gca,'box','on','FontName','times new roman','FontSize',12);
        
        set(gca,'color',[1 1 1]);
        
        set(gca,'xlim',[FPT_XMIN FPT_XMAX],'xtick',linspace(FPT_XMIN,FPT_XMAX,4));
        
        set(gca,'ylim',[FPT_YMIN FPT_YMAX],'ytick',linspace(FPT_YMIN,FPT_YMAX,4));
        
        set(gcf,'inverthardcopy','off'); %保存图片时按照设置的颜色而不自动调节
        
        set(gcf,'color',[1 1 1]);
        
        hold on 
     
        %---------绘制岸线、岛------------------------
        %         fill(land(:,1),land(:,2),bgcolor);
        plot(land(:,1),land(:,2),'color',[0.2 0.2 0.2]);
        hold on
        
        for k2 =1:nisland;
            %             eval(['fill( island' int2str(k2) '(:,1),island' int2str(k2) '(:,2),[' num2str(bgcolor) ']);']);
            eval(['plot( island' int2str(k2) '(:,1),island' int2str(k2) '(:,2),''color'',[0.2 0.2 0.2]);']);
            hold on
        end     
        %         fill(nanhuibiantan_weiken(:,1),nanhuibiantan_weiken(:,2),bgcolor);
        plot(nanhuibiantan_weiken(:,1),nanhuibiantan_weiken(:,2),'color',[0.2 0.2 0.2]);
        hold on
        
        plot(deep_water_way(:,1),deep_water_way(:,2),'linewidth',1,'color',[0.2 0.2 0.2]);
        hold on
        
        
        if (FPT_COORDINATE=='BL') %按照图像的坐标轴选取绘制相应的label
            xlabel('Longitude (\circE)');
            ylabel('Latitude (\circN)');
        elseif (FPT_COORDINATE=='XY')
            xlabel('Distance (km)');
            ylabel('Distance (km)');
        end     
        
        eval(['print(gcf,''-depsc'',''' fullfile(OPT_Pic_Path,'\Output') ''');']); %保存图像
        
    end
    
    
    
    
    
    %------------- Elevation Field Distribution Drawing -------------------
    if (LOG_EPT=='T') %如果需要画水位的平面图
        disp('Elevation Field Distribution Drawing ...')
        EPT_Pic_Path=fullfile(OUT_DIRE,'field_distri\elevation'); %水位平面图存放路径
        [s,mess,messid]=mkdir(EPT_Pic_Path); %建立水位平面图存放的文件夹
        EPT_Comp_Filepath=fullfile(IN_DIRE,'field_distri\elevation'); %模式的水位场计算结果文件路径
        
        
        EPT_TIME=str2num(EPT_TIME); %将读入的水位绘图时刻的字符数组转化为浮点数组
        Num_EPT_TIME=length(EPT_TIME); %总共需要画几个时刻的图像
        EPT_contourscale=linspace(-5,5,41); %需要画的等值线
        
        for i=1:Num_EPT_TIME %按照需要画图的数量进行循环
            close(figure(1));
            disp(['  The ' num2str(EPT_TIME(i)) ' Hours Drawing'])
            %------------------ Elevation Field Distribution Data Reading ------------------------
            EPT_Comp_Filename_SN=EPT_TIME(i)/(N_FPT/3600); %按照模式输出的时间间隔，计算需要画图时刻对应的序号
            EPT_Comp_Filename=sprintf('el_field_%06.6d',EPT_Comp_Filename_SN); %模式计算输出的水位文件
            
            if ~exist(fullfile(EPT_Comp_Filepath,EPT_Comp_Filename)) %判断水位场文件是否存在，不存在则输出提示信息
                disp(['    ' EPT_Comp_Filename ' dose not exist, drawing failed, please check the setting!']) %数据文件不存在，输出提示信息
                
            else %如果数据文件存在则进行绘图
                eval(['fid_EPT_Comp = fopen(''' fullfile(EPT_Comp_Filepath,EPT_Comp_Filename) ''',''r'',''b'');']); %打开文件，由于该二进制文件是用“BIG_ENDIAN”封装，所以打开是选择参数‘b’
                brecord=fread(fid_EPT_Comp,1,'integer*4'); %读记录大小信息
                thour=fread(fid_EPT_Comp,1,'real*4');
                EPT_Comp_Data=fread(fid_EPT_Comp,IM*JM,'real*4');
                erecord=fread(fid_EPT_Comp,1,'integer*4'); %读记录大小信息
                
                if (brecord~=erecord) %如果一条记录的前后所记录的大小不一致，则退出程序
                    disp(['Error in reading file ' fullfile(EPT_Comp_Filepath,EPT_Comp_Filename) ', please check!'])
                    break
                end
                EPT_Comp_Data=reshape(EPT_Comp_Data,IM,JM); %将水位数据转化为矩阵形式
                fclose(fid_EPT_Comp); %关闭计算输出的水位文件
                %---------------- End Elevation Field Distribution Data Reading ----------------------
                
                for i1=1:IM
                    for j1=1:JM
                        if (h(i1,j1)+EPT_Comp_Data(i1,j1)<DMIN||h(i1,j1)<-10) %将潮滩和陆地处的水位赋值为nan
                            EPT_Comp_Data(i1,j1)=NaN;
                        end
                    end
                end  
                
                if (FPT_COLOR_TYPE=='C') %根据绘图的颜色类型，选择颜色包文件
                    map=load ('Colormap\cm_el_C.dat');
                elseif (FPT_COLOR_TYPE=='G')
                    map=load ('Colormap\cm_el_G.dat');
                end
                map=[bgcolor;map]; %添加背景颜色，让水位为nan的位置按背景色画
                colormap(map);  %读入颜色包
                
                Pic_time=datevec(Init_day+thour/24); %根据模式设置的起始时间和该数据记录的thour计算实际的时间，Pic_time包含 年、月、日、时、分、秒 信息
                
                set(gca,'box','on','FontName','times new roman','FontSize',12);
                set(gca,'color',bgcolor);
                set(gca,'xlim',[FPT_XMIN FPT_XMAX],'xtick',linspace(FPT_XMIN,FPT_XMAX,4));
                set(gca,'ylim',[FPT_YMIN FPT_YMAX],'ytick',linspace(FPT_YMIN,FPT_YMAX,4));
                set(gcf,'inverthardcopy','off'); %保存图片时按照设置的颜色而不自动调节
                set(gcf,'color',[1 1 1]);
                hold on 
                
                warning off; %画contourf时会有很多警告信息，但不影响画图，所以关掉warning
                contourf(xr,yr,EPT_Comp_Data,EPT_contourscale);
                shading flat; 
                
                EPT_Comp_Data_Max=max(max(EPT_Comp_Data)); %求出水位数据的最大值最小值
                EPT_Comp_Data_Min=min(min(EPT_Comp_Data));
                cmax=ceil(max(abs(EPT_Comp_Data_Max),abs(EPT_Comp_Data_Min))); %求出最大的振幅，并向上取整
                cmin=-cmax; %因为颜色包里的颜色是对称的，所以要使正负坐标对称
                caxis([cmin-(cmax-cmin)/47 cmax]); %为了把背景色放到颜色包，而又不显示在colorbar里，‘47’是由于颜色包里有47个颜色
                
                if (FPT_COLOR_TYPE=='C') %如果画彩色图，则需要画colorbar
                    hc=colorbar;
                    set(hc,'ylim',[cmin cmax],'FontName','times new roman','FontSize',12); %设置colorbar的显示范围
                end
                
                [cc,hh]=contour(xr,yr,EPT_Comp_Data,EPT_contourscale);
                clabel(cc,hh,'FontSize',10,'FontName','Times new roman');
                set(hh,'EdgeColor',[0.3 0.3 0.3],'linewidth',1) %设置等值线的颜色和粗细
                
                % 				title([num2str(Pic_time(1)) '年' num2str(Pic_time(2)) '月' num2str(Pic_time(3)) '日' num2str(Pic_time(4)) '时 水位分布'])
                title([num2str(Pic_time(4)) ':00 ' mon(Pic_time(2),:) ' ' num2str(Pic_time(3)) ' ' num2str(Pic_time(1)) ' (GMT +8)'])    
                
                hold on  
                
                %---------绘制岸线、岛------------------------
                fill(land(:,1),land(:,2),bgcolor);
                plot(land(:,1),land(:,2),'color',[0.2 0.2 0.2]);
                
                hold on
                
                for k2 =1:nisland;
                    eval(['fill( island' int2str(k2) '(:,1),island' int2str(k2) '(:,2),[' num2str(bgcolor) ']);']);
                    eval(['plot( island' int2str(k2) '(:,1),island' int2str(k2) '(:,2),''color'',[0.2 0.2 0.2]);']);
                    hold on
                end         
                
                fill(nanhuibiantan_weiken(:,1),nanhuibiantan_weiken(:,2),bgcolor);
                plot(nanhuibiantan_weiken(:,1),nanhuibiantan_weiken(:,2),'color',[0.2 0.2 0.2]);
                hold on
                
                if (FPT_COORDINATE=='BL') %按照图像的坐标轴选取绘制相应的label
                    xlabel('Longitude (\circE)');
                    ylabel('Latitude (\circN)');
                elseif (FPT_COORDINATE=='XY')
                    xlabel('Distance (km)');
                    ylabel('Distance (km)');
                end
                
                eval(['print(gcf,''-dpng'',''' fullfile(EPT_Pic_Path,EPT_Comp_Filename) ''');']); %保存图像
                
            end
            
        end      
    end    
    %------------- End Elevation Field Distribution Drawing --------------
    
    %------------- Velocity Field Distribution Drawing -------------------
    if (LOG_VPT_UV=='T'||LOG_VPT_SD=='T') %如果需要画流场的平面图
        disp('Velocity Field Distribution Drawing ...')
        VPT_Pic_Path=fullfile(OUT_DIRE,'field_distri\current'); %盐度平面图存放路径
        [s,mess,messid]=mkdir(VPT_Pic_Path); %建立盐度平面图的存放文件夹
        EPT_Comp_Filepath=fullfile(IN_DIRE,'field_distri\elevation'); %模式的水位场计算结果文件路径
        VPT_Comp_Filepath=fullfile(IN_DIRE,'field_distri\current'); %模式的盐度场计算结果文件路径
        
        
        VPT_TIME=str2num(VPT_TIME); %将读入的盐度绘图时刻的字符数组转化为浮点数组
        Num_VPT_TIME=length(VPT_TIME); %总共需要画几个时刻的图像
        VPT_LAYER=str2num(VPT_LAYER); %将读入的盐度绘图层数的字符数组转化为浮点数组
        Num_VPT_LAYER=length(VPT_LAYER); %总共需要画几个层次的图像
        H_contourscale=[-1 0 5 10 15 20 25 30 50 70 100]; %矢量彩图中水深等值线设置
        SD_contourscale=[-1 0 5 10 15 20 25 30 50 70 100]; %流速量值等值线设置
        
        
        
        for i=1:Num_VPT_TIME %按照需要画图的数量进行循环
            
            disp(['  The ' num2str(VPT_TIME(i)) ' Hours Drawing'])
            
            %------ Elevation Field Distribution Data Reading ,For Check Tideflat Grid ----------
            EPT_Comp_Filename_SN=VPT_TIME(i)/(N_FPT/3600); %按照模式输出的时间间隔，计算需要画图时刻对应的序号
            EPT_Comp_Filename=sprintf('el_field_%06.6d',EPT_Comp_Filename_SN); %模式计算输出的水位文件
            VPT_Comp_Filename_SN=VPT_TIME(i)/(N_FPT/3600); %按照模式输出的时间间隔，计算需要画图时刻对应的序号
            VPT_Comp_Filename=sprintf('v_field_%06.6d',VPT_Comp_Filename_SN); %模式计算输出的水位文件           
            
            if ~(exist(fullfile(EPT_Comp_Filepath,EPT_Comp_Filename))&&exist(fullfile(VPT_Comp_Filepath,VPT_Comp_Filename))) %判断水位场文件和流速场文件是否都存在，不存在则输出提示信息
                disp(['    ' EPT_Comp_Filename ' or ' VPT_Comp_Filename ' dose not exist, drawing failed, please check the setting!']) %数据文件不存在，输出提示信息
            else %数据文件存在，则绘图
                eval(['fid_EPT_Comp = fopen(''' fullfile(EPT_Comp_Filepath,EPT_Comp_Filename) ''',''r'',''b'');']); %打开文件，由于该二进制文件是用“BIG_ENDIAN”封装，所以打开是选择参数‘b’
                brecord=fread(fid_EPT_Comp,1,'integer*4'); %读记录大小信息
                thour=fread(fid_EPT_Comp,1,'real*4');
                EPT_Comp_Data=fread(fid_EPT_Comp,IM*JM,'real*4');
                erecord=fread(fid_EPT_Comp,1,'integer*4'); %读记录大小信息
                
                if (brecord~=erecord) %如果一条记录的前后所记录的大小不一致，则退出程序
                    disp(['Error in reading file ' fullfile(EPT_Comp_Filepath,EPT_Comp_Filename) ', please check!'])
                    break
                end
                EPT_Comp_Data=reshape(EPT_Comp_Data,IM,JM); %将水位数据转化为矩阵形式
                fclose(fid_EPT_Comp); %关闭计算输出的水位文件
                %---------------- End Elevation Field Distribution Data Reading ----------------------
                
                
                
                %--------------------- Velocity Field Distribution Data Reading -----------------------
                eval(['fid_VPT_Comp = fopen(''' fullfile(VPT_Comp_Filepath,VPT_Comp_Filename) ''',''r'',''b'');']); %打开文件，由于该二进制文件是用“BIG_ENDIAN”封装，所以打开是选择参数‘b’
                brecord=fread(fid_VPT_Comp,1,'integer*4'); %读记录大小信息
                thour=fread(fid_VPT_Comp,1,'real*4');
                VPT_Comp_U=fread(fid_VPT_Comp,IM*JM*KB,'real*4');
                VPT_Comp_V=fread(fid_VPT_Comp,IM*JM*KB,'real*4');
                erecord=fread(fid_VPT_Comp,1,'integer*4'); %读记录大小信息
                
                if (brecord~=erecord) %如果一条记录的前后所记录的大小不一致，则退出程序
                    disp(['Error in reading file ' fullfile(VPT_Comp_Filepath,VPT_Comp_Filename) ', please check!'])
                    break
                end
                VPT_Comp_U=reshape(VPT_Comp_U,IM,JM,KB); %将流速U转化为矩阵形式
                VPT_Comp_V=reshape(VPT_Comp_V,IM,JM,KB); %将流速V转化为矩阵形式
                fclose(fid_VPT_Comp); %关闭计算输出的流速文件
                %------------------- End Velocity Field Distribution Data Reading --------------------
                
                
                
                %初始化用于绘制流速矢量图的变量，都赋值为nan
                X_VPT=NaN*ones(IM,JM);
                Y_VPT=NaN*ones(IM,JM);
                U_VPT=NaN*ones(IM,JM,KB);
                V_VPT=NaN*ones(IM,JM,KB);
                
                %按照设置的绘图间隔VPT_INTERVAL，将有效点的位置坐标和流速赋值
                for i2=1:VPT_INTERVAL:IM 
                    for j2=1:VPT_INTERVAL:JM
                        X_VPT(i2,j2)=xr(i2,j2);
                        Y_VPT(i2,j2)=yr(i2,j2);                       
                        U_VPT(i2,j2,:)=VPT_Comp_U(i2,j2,:);
                        V_VPT(i2,j2,:)=VPT_Comp_V(i2,j2,:);
                    end
                end    				
                
                
                
                H_VPT=h; %将水深数据赋值给一个临时变量H_VPT，可以将H_VPT进行操作
                for i1=1:IM
                    for j1=1:JM
                        if (h(i1,j1)+EPT_Comp_Data(i1,j1)<DMIN||h(i1,j1)<-10) %将潮滩和陆地处的水深和流速赋值为nan
                            H_VPT(i1,j1)=NaN;
                            X_VPT(i1,j1)=NaN;
                            Y_VPT(i1,j1)=NaN;                    
                            U_VPT(i1,j1,:)=NaN; %用于画矢量图
                            V_VPT(i1,j1,:)=NaN;  
                            VPT_Comp_U(i1,j1,:)=NaN; %用于画量值图
                            VPT_Comp_V(i1,j1,:)=NaN;  
                            
                        end
                    end
                end  
                
                %为了画出流矢大小标记，在数据末尾添加一个水平方向的数据
                X_VPT(i1,j1+1)=FPT_XMIN+(FPT_XMAX-FPT_XMIN)*0.05; %横向为从左向右5%的位置
                Y_VPT(i1,j1+1)=FPT_YMIN+(FPT_YMAX-FPT_YMIN)*0.1; %横向为从小而上10%的位置
                U_VPT(i1,j1+1,:)=VPT_SCALE;
                V_VPT(i1,j1+1,:)=0;
                
                %为了画出流矢大小标记，在数据末尾添加一个垂直方向的数据
                X_VPT(i1+1,j1+1)=FPT_XMIN+(FPT_XMAX-FPT_XMIN)*0.05; %横向为从左向右5%的位置
                Y_VPT(i1+1,j1+1)=FPT_YMIN+(FPT_YMAX-FPT_YMIN)*0.12; %横向为从小而上12%的位置
                U_VPT(i1+1,j1+1,:)=0;
                V_VPT(i1+1,j1+1,:)=VPT_SCALE;		
                
                %在quiver函数中设定箭头按实际大小画，而不是自动调节大小，所以需要把将原始数据放大到和坐标轴相应的大小尺度
                U_VPT=U_VPT*(FPT_XMAX-FPT_XMIN)*0.05/VPT_SCALE; %(FPT_XMAX-FPT_XMIN)*0.05/VPT_SCALE表示将x轴的5%的长度对应VPT_SCALE的大小
                V_VPT=V_VPT*(FPT_YMAX-FPT_YMIN)*0.05/VPT_SCALE; %(FPT_YMAX-FPT_YMIN)*0.05/VPT_SCALE表示将y轴的5%的长度对应VPT_SCALE的大小
                
                Pic_time=datevec(Init_day+thour/24); %根据模式设置的起始时间和该数据记录的thour计算实际的时间，Pic_time包含 年、月、日、时、分、秒 信息
                
                
                for k=1:Num_VPT_LAYER
                    
                    if (LOG_VPT_UV=='T') %绘制流速平面矢量图
                        
                        close(figure(1));
                        disp(['    Layer ' num2str(VPT_LAYER(k))])
                        if (FPT_COLOR_TYPE=='C') %根据绘图的颜色类型，选择颜色包文件
                            map=load ('Colormap\cm_v_C.dat');
                        elseif (FPT_COLOR_TYPE=='G')
                            map=load ('Colormap\cm_v_G.dat');
                        end
                        map=[bgcolor;map]; %添加背景颜色，让盐度为nan的位置按背景色画
                        colormap(map);  %读入颜色包
                        
                        
                        set(gca,'box','on','FontName','times new roman','FontSize',12);
                        set(gca,'color',bgcolor);
                        set(gca,'xlim',[FPT_XMIN FPT_XMAX],'xtick',linspace(FPT_XMIN,FPT_XMAX,4));
                        set(gca,'ylim',[FPT_YMIN FPT_YMAX],'ytick',linspace(FPT_YMIN,FPT_YMAX,4));
                        set(gcf,'inverthardcopy','off'); %保存图片时按照设置的颜色而不自动调节
                        set(gcf,'color',[1 1 1]);
                        hold on 
                        
                        
                        warning off; %画contourf时会有很多警告信息，但不影响画图，所以关掉warning
                        contourf(xr,yr,H_VPT,H_contourscale);
                        shading flat; 
                        caxis([-1 70]); 
                        if (FPT_COLOR_TYPE=='C') %如果画彩色图，则需要画colorbar
                            hc=colorbar;
                            set(hc,'ylim',[0 70],'FontName','times new roman','FontSize',12); %设置colorbar的显示范围
                            [cc,hh]=contour(xr,yr,H_VPT,H_contourscale); %彩图画水深的等值线，灰度图不画
                            clabel(cc,hh,'FontSize',10,'FontName','Times new roman');
                            set(hh,'EdgeColor',[0.3 0.3 0.3],'linewidth',1) %设置等值线的颜色和粗细
                        end
                        
                        eval(['quiver(X_VPT,Y_VPT,U_VPT(:,:,' num2str(VPT_LAYER(k)) '),V_VPT(:,:,' num2str(VPT_LAYER(k)) '),0,''k'');']);
                        
                        % 					if (VPT_LAYER(k)==1)
                        % 						VPT_LAYER_NAME='表';
                        % 					elseif (VPT_LAYER(k)==5)
                        % 						VPT_LAYER_NAME='底';
                        % 					else
                        % 						VPT_LAYER_NAME=['第' num2str(VPT_LAYER(k))];
                        % 					end
                        % 					title([num2str(Pic_time(1)) '年' num2str(Pic_time(2)) '月' num2str(Pic_time(3)) '日' num2str(Pic_time(4)) '时 ' VPT_LAYER_NAME '层流速分布'])
                        
                        if (VPT_LAYER(k)==1)
                            VPT_LAYER_NAME='Surface';
                        elseif (VPT_LAYER(k)==KB-1)
                            VPT_LAYER_NAME='Bottom';
                        elseif (VPT_LAYER(k)==2)
                            VPT_LAYER_NAME='2nd layer';
                        elseif (VPT_LAYER(k)==3)
                            VPT_LAYER_NAME='3rd layer';
                        else
                            VPT_LAYER_NAME=[num2str(VPT_LAYER(k)) 'th layer'];
                        end
                        title([num2str(Pic_time(4)) ':00 ' mon(Pic_time(2),:) ' ' num2str(Pic_time(3)) ' ' num2str(Pic_time(1)) ' (GMT +8)   ' VPT_LAYER_NAME])
                        
                        hold on  
                        
                        %---------绘制岸线、岛------------------------
                        fill(land(:,1),land(:,2),bgcolor);
                        plot(land(:,1),land(:,2),'color',[0.2 0.2 0.2]);
                        
                        hold on
                        
                        for k2 =1:nisland;
                            eval(['fill( island' int2str(k2) '(:,1),island' int2str(k2) '(:,2),[' num2str(bgcolor) ']);']);
                            eval(['plot( island' int2str(k2) '(:,1),island' int2str(k2) '(:,2),''color'',[0.2 0.2 0.2]);']);
                            hold on
                        end         
                        
                        fill(nanhuibiantan_weiken(:,1),nanhuibiantan_weiken(:,2),bgcolor);
                        plot(nanhuibiantan_weiken(:,1),nanhuibiantan_weiken(:,2),'color',[0.2 0.2 0.2]);
                        hold on
                        
                        text(FPT_XMIN+(FPT_XMAX-FPT_XMIN)*0.05,FPT_YMIN+(FPT_YMAX-FPT_YMIN)*0.075,[num2str(VPT_SCALE) 'm/s'],'FontName','times new roman','FontSize',10);
                        
                        if (FPT_COORDINATE=='BL') %按照图像的坐标轴选取绘制相应的label
                            xlabel('Longitude (\circE)');
                            ylabel('Latitude (\circN)');
                        elseif (FPT_COORDINATE=='XY')
                            xlabel('Distance (km)');
                            ylabel('Distance (km)');
                        end
                        
                        eval(['print(gcf,''-dpng'',''' fullfile(VPT_Pic_Path,VPT_Comp_Filename) '_UV_' num2str(VPT_LAYER(k)) ''');']); %保存图像
                        
                    end %绘制流速平面矢量图 end
                    
                    
                    
                    
                    if (LOG_VPT_SD=='T') %绘制流速平面量值图（等值线图）
                        
                        close(figure(1));
                        disp(['    Layer ' num2str(VPT_LAYER(k))])
                        if (FPT_COLOR_TYPE=='C') %根据绘图的颜色类型，选择颜色包文件
                            map=load ('Colormap\cm_speed_C.dat');
                        elseif (FPT_COLOR_TYPE=='G')
                            map=load ('Colormap\cm_speed_G.dat');
                        end
                        map=[bgcolor;map]; %添加背景颜色，让盐度为nan的位置按背景色画
                        colormap(map);  %读入颜色包
                        
                        
                        set(gca,'box','on','FontName','times new roman','FontSize',12);
                        set(gca,'color',bgcolor);
                        set(gca,'xlim',[FPT_XMIN FPT_XMAX],'xtick',linspace(FPT_XMIN,FPT_XMAX,4));
                        set(gca,'ylim',[FPT_YMIN FPT_YMAX],'ytick',linspace(FPT_YMIN,FPT_YMAX,4));
                        set(gcf,'inverthardcopy','off'); %保存图片时按照设置的颜色而不自动调节
                        set(gcf,'color',[1 1 1]);
                        hold on 
                        
                        VPT_Comp_Data(:,:,VPT_LAYER(k))=sqrt(VPT_Comp_U(:,:,VPT_LAYER(k)).^2+VPT_Comp_V(:,:,VPT_LAYER(k)).^2); %把uv转化为流速
                        
                        warning off; %画contourf时会有很多警告信息，但不影响画图，所以关掉warning
                        contourf(xr,yr,VPT_Comp_Data(:,:,VPT_LAYER(k)));
                        
                        shading flat; 
                        
                        
                        % 						caxis([-6/64 3]); %给定最大量值为3，或是按照整个计算区域最大值指定最大值
                        VPT_Comp_Data_Max=max(max(VPT_Comp_Data(:,:,VPT_LAYER(k)))); %求出流速数据的最大值最小值
                        VPT_Comp_Data_Min=min(min(VPT_Comp_Data(:,:,VPT_LAYER(k))));
                        cmax=ceil(max(abs(VPT_Comp_Data_Max),abs(VPT_Comp_Data_Min))); %求出最大的振幅，并向上取整
                        caxis([-cmax/64 cmax]); %为了把背景色放到颜色包，而又不显示在colorbar里，‘64’是由于颜色包里有64个颜色                 
                        
                        if (FPT_COLOR_TYPE=='C') %如果画彩色图，则需要画colorbar
                            hc=colorbar;
                            set(hc,'ylim',[0 cmax],'FontName','times new roman','FontSize',12); %设置colorbar的显示范围
                            [cc,hh]=contour(xr,yr,VPT_Comp_Data(:,:,VPT_LAYER(k))); 
                            clabel(cc,hh,'FontSize',10,'FontName','Times new roman');
                            set(hh,'EdgeColor',[0.3 0.3 0.3],'linewidth',1) %设置等值线的颜色和粗细
                        elseif (FPT_COLOR_TYPE=='G')
                            [cc,hh]=contour(xr,yr,VPT_Comp_Data(:,:,VPT_LAYER(k))); 
                            clabel(cc,hh,'FontSize',10,'FontName','Times new roman');
                            set(hh,'EdgeColor',[0.3 0.3 0.3],'linewidth',1) %设置等值线的颜色和粗细
                        end
                        
                        % 						eval(['quiver(X_VPT,Y_VPT,U_VPT(:,:,' num2str(VPT_LAYER(k)) '),V_VPT(:,:,' num2str(VPT_LAYER(k)) '),0,''k'');']);
                        
                        % 					if (VPT_LAYER(k)==1)
                        % 						VPT_LAYER_NAME='表';
                        % 					elseif (VPT_LAYER(k)==5)
                        % 						VPT_LAYER_NAME='底';
                        % 					else
                        % 						VPT_LAYER_NAME=['第' num2str(VPT_LAYER(k))];
                        % 					end
                        % 					title([num2str(Pic_time(1)) '年' num2str(Pic_time(2)) '月' num2str(Pic_time(3)) '日' num2str(Pic_time(4)) '时 ' VPT_LAYER_NAME '层流速分布'])
                        
                        if (VPT_LAYER(k)==1)
                            VPT_LAYER_NAME='Surface';
                        elseif (VPT_LAYER(k)==KB-1)
                            VPT_LAYER_NAME='Bottom';
                        elseif (VPT_LAYER(k)==2)
                            VPT_LAYER_NAME='2nd layer';
                        elseif (VPT_LAYER(k)==3)
                            VPT_LAYER_NAME='3rd layer';
                        else
                            VPT_LAYER_NAME=[num2str(VPT_LAYER(k)) 'th layer'];
                        end
                        title([num2str(Pic_time(4)) ':00 ' mon(Pic_time(2),:) ' ' num2str(Pic_time(3)) ' ' num2str(Pic_time(1)) ' (GMT +8)   ' VPT_LAYER_NAME])
                        
                        hold on  
                        
                        %---------绘制岸线、岛------------------------
                        fill(land(:,1),land(:,2),bgcolor);
                        plot(land(:,1),land(:,2),'color',[0.2 0.2 0.2]);
                        
                        hold on
                        
                        for k2 =1:nisland;
                            eval(['fill( island' int2str(k2) '(:,1),island' int2str(k2) '(:,2),[' num2str(bgcolor) ']);']);
                            eval(['plot( island' int2str(k2) '(:,1),island' int2str(k2) '(:,2),''color'',[0.2 0.2 0.2]);']);
                            hold on
                        end         
                        
                        fill(nanhuibiantan_weiken(:,1),nanhuibiantan_weiken(:,2),bgcolor);
                        plot(nanhuibiantan_weiken(:,1),nanhuibiantan_weiken(:,2),'color',[0.2 0.2 0.2]);
                        hold on
                        
                        % 						text(FPT_XMIN+(FPT_XMAX-FPT_XMIN)*0.05,FPT_YMIN+(FPT_YMAX-FPT_YMIN)*0.075,[num2str(VPT_SCALE) 'm/s'],'FontName','times new roman','FontSize',10);
                        
                        if (FPT_COORDINATE=='BL') %按照图像的坐标轴选取绘制相应的label
                            xlabel('Longitude (\circE)');
                            ylabel('Latitude (\circN)');
                        elseif (FPT_COORDINATE=='XY')
                            xlabel('Distance (km)');
                            ylabel('Distance (km)');
                        end
                        
                        eval(['print(gcf,''-dpng'',''' fullfile(VPT_Pic_Path,VPT_Comp_Filename) '_SD_' num2str(VPT_LAYER(k)) ''');']); %保存图像
                        
                    end %绘制流速平面量值图（等值线图）end
                    
                end
                clear H_VPT X_VPT Y_VPT U_VPT V_VPT
            end
        end
    end
    %------------- End Velocity Field Distribution Drawing ----------------    
    
    %------------- Salinity Field Distribution Drawing -------------------
    if (LOG_SPT=='T') %如果需要画盐度的平面图
        disp('Salinity Field Distribution Drawing ...')
        SPT_Pic_Path=fullfile(OUT_DIRE,'field_distri\salinity'); %盐度平面图存放路径
        [s,mess,messid]=mkdir(SPT_Pic_Path); %建立盐度平面图的存放文件夹
        EPT_Comp_Filepath=fullfile(IN_DIRE,'field_distri\elevation'); %模式的水位场计算结果文件路径
        SPT_Comp_Filepath=fullfile(IN_DIRE,'field_distri\salinity'); %模式的盐度场计算结果文件路径
        
        
        SPT_TIME=str2num(SPT_TIME); %将读入的盐度绘图时刻的字符数组转化为浮点数组
        Num_SPT_TIME=length(SPT_TIME); %总共需要画几个时刻的图像
        SPT_LAYER=str2num(SPT_LAYER); %将读入的盐度绘图层数的字符数组转化为浮点数组
        Num_SPT_LAYER=length(SPT_LAYER); %总共需要画几个层次的图像
        SPT_contourscale=[-0.5 0 0.5 1 3 5 10 15 20 25 30 35]; %需要画的等值线
        
        
        
        
        for i=1:Num_SPT_TIME %按照需要画图的数量进行循环
            
            disp(['  The ' num2str(SPT_TIME(i)) ' Hours Drawing'])
            
            %------ Elevation Field Distribution Data Reading ,For Check Tideflat Grid ----------
            EPT_Comp_Filename_SN=SPT_TIME(i)/(N_FPT/3600); %按照模式输出的时间间隔，计算需要画图时刻对应的序号
            EPT_Comp_Filename=sprintf('el_field_%06.6d',EPT_Comp_Filename_SN); %模式计算输出的水位文件
            SPT_Comp_Filename_SN=SPT_TIME(i)/(N_FPT/3600); %按照模式输出的时间间隔，计算需要画图时刻对应的序号
            SPT_Comp_Filename=sprintf('s_field_%06.6d',SPT_Comp_Filename_SN); %模式计算输出的水位文件           
            
            if ~(exist(fullfile(EPT_Comp_Filepath,EPT_Comp_Filename))&&exist(fullfile(SPT_Comp_Filepath,SPT_Comp_Filename))) %判断水位场文件和盐度场文件是否都存在，不存在则输出提示信息
                disp(['    ' EPT_Comp_Filename ' or ' SPT_Comp_Filename ' dose not exist, drawing failed, please check the setting!']) %数据文件不存在，输出提示信息
            else %数据文件存在，则绘图
                eval(['fid_EPT_Comp = fopen(''' fullfile(EPT_Comp_Filepath,EPT_Comp_Filename) ''',''r'',''b'');']); %打开文件，由于该二进制文件是用“BIG_ENDIAN”封装，所以打开是选择参数‘b’
                brecord=fread(fid_EPT_Comp,1,'integer*4'); %读记录大小信息
                thour=fread(fid_EPT_Comp,1,'real*4');
                EPT_Comp_Data=fread(fid_EPT_Comp,IM*JM,'real*4');
                erecord=fread(fid_EPT_Comp,1,'integer*4'); %读记录大小信息
                
                if (brecord~=erecord) %如果一条记录的前后所记录的大小不一致，则退出程序
                    disp(['Error in reading file ' fullfile(EPT_Comp_Filepath,EPT_Comp_Filename) ', please check!'])
                    break
                end
                EPT_Comp_Data=reshape(EPT_Comp_Data,IM,JM); %将水位数据转化为矩阵形式
                fclose(fid_EPT_Comp); %关闭计算输出的水位文件
                %---------------- End Elevation Field Distribution Data Reading ----------------------
                
                
                
                %--------------------- Salinity Field Distribution Data Reading -----------------------
                eval(['fid_SPT_Comp = fopen(''' fullfile(SPT_Comp_Filepath,SPT_Comp_Filename) ''',''r'',''b'');']); %打开文件，由于该二进制文件是用“BIG_ENDIAN”封装，所以打开是选择参数‘b’
                brecord=fread(fid_SPT_Comp,1,'integer*4'); %读记录大小信息
                thour=fread(fid_SPT_Comp,1,'real*4');
                SPT_Comp_Data=fread(fid_SPT_Comp,IM*JM*KB,'real*4');
                erecord=fread(fid_SPT_Comp,1,'integer*4'); %读记录大小信息
                
                if (brecord~=erecord) %如果一条记录的前后所记录的大小不一致，则退出程序
                    disp(['Error in reading file ' fullfile(SPT_Comp_Filepath,SPT_Comp_Filename) ', please check!'])
                    break
                end
                SPT_Comp_Data=reshape(SPT_Comp_Data,IM,JM,KB); %将盐度转化为矩阵形式
                fclose(fid_SPT_Comp); %关闭计算输出的水位文件
                %------------------- End Salinity Field Distribution Data Reading -----------------------
                
                
                for i1=1:IM
                    for j1=1:JM
                        if (h(i1,j1)+EPT_Comp_Data(i1,j1)<DMIN||h(i1,j1)<-10) %将潮滩和陆地处的盐度赋值为nan
                            SPT_Comp_Data(i1,j1,:)=NaN;
                        end
                        
                        for k=1:KB-1
                            if SPT_Comp_Data(i1,j1,k)<0
                                SPT_Comp_Data(i1,j1,k)=0;
                            end
                        end
                    end
                end  
                
                
                
                Pic_time=datevec(Init_day+thour/24); %根据模式设置的起始时间和该数据记录的thour计算实际的时间，Pic_time包含 年、月、日、时、分、秒 信息
                
                
                for k=1:Num_SPT_LAYER
                    close(figure(1));
                    disp(['    Layer ' num2str(SPT_LAYER(k))])
                    if (FPT_COLOR_TYPE=='C') %根据绘图的颜色类型，选择颜色包文件
                        map=load ('Colormap\cm_s_C.dat');
                    elseif (FPT_COLOR_TYPE=='G')
                        map=load ('Colormap\cm_s_G.dat');
                    end
                    map=[bgcolor;map]; %添加背景颜色，让盐度为nan的位置按背景色画
                    colormap(map);  %读入颜色包
                    
                    
                    set(gca,'box','on','FontName','times new roman','FontSize',12);
                    set(gca,'color',bgcolor);
                    set(gca,'xlim',[FPT_XMIN FPT_XMAX],'xtick',linspace(FPT_XMIN,FPT_XMAX,4));
                    set(gca,'ylim',[FPT_YMIN FPT_YMAX],'ytick',linspace(FPT_YMIN,FPT_YMAX,4));
                    set(gcf,'inverthardcopy','off'); %保存图片时按照设置的颜色而不自动调节
                    set(gcf,'color',[1 1 1]);
                    hold on 
                    
                    
                    warning off; %画contourf时会有很多警告信息，但不影响画图，所以关掉warning
                    eval(['contourf(xr,yr,SPT_Comp_Data(:,:,' num2str(SPT_LAYER(k)) '),SPT_contourscale);']);
                    shading flat; 
                    caxis([-0.5 35]); %为了把背景色放到颜色包，而又不显示在colorbar里，‘/20’是由于颜色包里有20个颜色
                    if (FPT_COLOR_TYPE=='C') %如果画彩色图，则需要画colorbar
                        hc=colorbar;
                        set(hc,'ylim',[0 35],'FontName','times new roman','FontSize',12); %设置colorbar的显示范围
                    end
                    
                    
                    eval(['[cc,hh]=contour(xr,yr,SPT_Comp_Data(:,:,' num2str(SPT_LAYER(k)) '),SPT_contourscale);']);
                    clabel(cc,hh,'FontSize',10,'FontName','Times new roman');
                    set(hh,'EdgeColor',[0.3 0.3 0.3],'linewidth',1) %设置等值线的颜色和粗细
 
                    if (SPT_LAYER(k)==1)
                        SPT_LAYER_NAME='Surface';
                    elseif (SPT_LAYER(k)==KB-1)
                        SPT_LAYER_NAME='Bottom';
                    elseif (SPT_LAYER(k)==2)
                        SPT_LAYER_NAME='2nd layer';
                    elseif (SPT_LAYER(k)==3)
                        SPT_LAYER_NAME='3rd layer';
                    else
                        SPT_LAYER_NAME=[num2str(SPT_LAYER(k)) 'th layer'];
                    end
                    title([num2str(Pic_time(4)) ':00 ' mon(Pic_time(2),:) ' ' num2str(Pic_time(3)) ' ' num2str(Pic_time(1)) ' (GMT +8)   ' SPT_LAYER_NAME])
                    hold on  
                    
                    %---------绘制岸线、岛------------------------
                    fill(land(:,1),land(:,2),bgcolor);
                    plot(land(:,1),land(:,2),'color',[0.2 0.2 0.2]);
                    
                    hold on
                    
                    for k2 =1:nisland;
                        eval(['fill( island' int2str(k2) '(:,1),island' int2str(k2) '(:,2),[' num2str(bgcolor) ']);']);
                        eval(['plot( island' int2str(k2) '(:,1),island' int2str(k2) '(:,2),''color'',[0.2 0.2 0.2]);']);
                        hold on
                    end         
                    
                    fill(nanhuibiantan_weiken(:,1),nanhuibiantan_weiken(:,2),bgcolor);
                    plot(nanhuibiantan_weiken(:,1),nanhuibiantan_weiken(:,2),'color',[0.2 0.2 0.2]);
                    hold on
                    
                    if (FPT_COORDINATE=='BL') %按照图像的坐标轴选取绘制相应的label
                        xlabel('Longitude (\circE)');
                        ylabel('Latitude (\circN)');
                    elseif (FPT_COORDINATE=='XY')
                        xlabel('Distance (km)');
                        ylabel('Distance (km)');
                    end
                    
                    eval(['print(gcf,''-dpng'',''' fullfile(SPT_Pic_Path,SPT_Comp_Filename) '_' num2str(SPT_LAYER(k)) ''');']); %保存图像
                    
                end
            end
        end      
    end 
    %------------- End Salinity Field Distribution Drawing ----------------
    
    
    %----------- Residual Salinity Field Distribution Drawing -------------
    if (LOG_RSPT=='T') %如果需要画盐度的平面图
        disp('Residual Salinity Field Distribution Drawing ...')
        RSPT_Pic_Path=fullfile(OUT_DIRE,'residual_distri\field\salinity'); %余盐度平面图存放路径
        [s,mess,messid]=mkdir(RSPT_Pic_Path); %建立盐度平面图的存放文件夹
        RSPT_Comp_Filepath=fullfile(IN_DIRE,'resi_flux'); %模式的余盐度场计算结果文件路径
        
        RSPT_LAYER=str2num(RSPT_LAYER); %将读入的盐度绘图层数的字符数组转化为浮点数组
        Num_RSPT_LAYER=length(RSPT_LAYER); %总共需要画几个层次的图像
        RSPT_contourscale=[-0.5 0 0.5 1 3 5 10 15 20 25 30 35]; %需要画的等值线
        
        RSPT_Comp_Fileinfo=dir(fullfile(RSPT_Comp_Filepath,'resi_salt_flux_3d_*.out')); %模式计算结果文件信息（包括文件名，修改时间，大小，是否为文件路径）
        
        
        for i=1:length(RSPT_Comp_Fileinfo) %按照需要画图的数量进行循环
            
            %--------------------- Residual Salinity Field Distribution Data Reading -----------------------
            RSPT_Comp_Filename=RSPT_Comp_Fileinfo(i).name; %模式计算输出的余盐度文件           
            eval(['fid_RSPT_Comp = fopen(''' fullfile(RSPT_Comp_Filepath,RSPT_Comp_Filename) ''',''r'');']); %打开余盐度文件
            RSPT_Comp_Data=nan*ones(IM,JM,KB);
            Note=fscanf(fid_RSPT_Comp,'%s',7);
            RSPT_B=fscanf(fid_RSPT_Comp,'%f',1);
            Note=fscanf(fid_RSPT_Comp,'%s',2);
            RSPT_E=fscanf(fid_RSPT_Comp,'%f',1);
            Note=fgetl(fid_RSPT_Comp);
            Note=fgetl(fid_RSPT_Comp);
            disp(['  Residual Salinty Distribution from ' num2str(RSPT_B) ' to ' num2str(RSPT_E) ' Hour Drawing'])
            
            while(~feof(fid_RSPT_Comp)) %判读是否读到文件尾部
                I=fscanf(fid_RSPT_Comp,'%d',1); %读I
                J=fscanf(fid_RSPT_Comp,'%d',1); %读J
                Note=fscanf(fid_RSPT_Comp,'%f',2); %读X,Y
                RSPT_Comp_Data(I,J,1:KB-1)=fscanf(fid_RSPT_Comp,'%f',5); %读X,Y位置kbm1层的盐度值
                Note=fgetl(fid_RSPT_Comp); %读掉末尾的回车
            end
            
            fclose(fid_RSPT_Comp); %关闭计算输出的余盐度文件
            %------------------- End Residual Salinity Field Distribution Data Reading -----------------------
            
            for i1=1:IM
                for j1=1:JM
                    if (h(i1,j1)<RS_MIN_DEP) %小于RS_MIN_DEP位置不画
                        RSPT_Comp_Data(i1,j1,:)=NaN; 
                        
                        for k=1:KB-1
                            if RSPT_Comp_Data(i1,j1,k)<0
                                RSPT_Comp_Data(i1,j1,k)=0;
                            end
                        end
                    end
                end
            end  
            
            
            
            
            Pic_time_b=datevec(Init_day+RSPT_B/24); %根据模式设置的起始时间计算余盐度场统计的起始时间
            Pic_time_e=datevec(Init_day+RSPT_E/24); %根据模式设置的起始时间计算余盐度场统计的结束时间
            
            
            for k=1:Num_RSPT_LAYER
                close(figure(1));
                disp(['    Layer ' num2str(RSPT_LAYER(k))])
                if (FPT_COLOR_TYPE=='C') %根据绘图的颜色类型，选择颜色包文件
                    map=load ('Colormap\cm_s_C.dat');
                elseif (FPT_COLOR_TYPE=='G')
                    map=load ('Colormap\cm_s_G.dat');
                end
                map=[bgcolor;map]; %添加背景颜色，让盐度为nan的位置按背景色画
                colormap(map);  %读入颜色包
                
                
                set(gca,'box','on','FontName','times new roman','FontSize',12);
                set(gca,'color',bgcolor);
                set(gca,'xlim',[FPT_XMIN FPT_XMAX],'xtick',linspace(FPT_XMIN,FPT_XMAX,4));
                set(gca,'ylim',[FPT_YMIN FPT_YMAX],'ytick',linspace(FPT_YMIN,FPT_YMAX,4));
                set(gcf,'inverthardcopy','off'); %保存图片时按照设置的颜色而不自动调节
                set(gcf,'color',[1 1 1]);
                hold on 
                
                
                warning off; %画contourf时会有很多警告信息，但不影响画图，所以关掉warning
                eval(['contourf(xr,yr,RSPT_Comp_Data(:,:,' num2str(RSPT_LAYER(k)) '),RSPT_contourscale);']);
                shading flat; 
                caxis([-0.5 35]); %为了把背景色放到颜色包，而又不显示在colorbar里，‘/20’是由于颜色包里有20个颜色
                if (FPT_COLOR_TYPE=='C') %如果画彩色图，则需要画colorbar
                    hc=colorbar;
                    set(hc,'ylim',[0 35],'FontName','times new roman','FontSize',12); %设置colorbar的显示范围
                end
                
                
                eval(['[cc,hh]=contour(xr,yr,RSPT_Comp_Data(:,:,' num2str(RSPT_LAYER(k)) '),RSPT_contourscale);']);
                clabel(cc,hh,'FontSize',10,'FontName','Times new roman');
                set(hh,'EdgeColor',[0.3 0.3 0.3],'linewidth',1) %设置等值线的颜色和粗细
                
                
                
                if (RSPT_LAYER(k)==1)
                    RSPT_LAYER_NAME='Surface';
                elseif (RSPT_LAYER(k)==KB-1)
                    RSPT_LAYER_NAME='Bottom';
                elseif (RSPT_LAYER(k)==2)
                    RSPT_LAYER_NAME='2nd layer';
                elseif (RSPT_LAYER(k)==3)
                    RSPT_LAYER_NAME='3rd layer';
                else
                    RSPT_LAYER_NAME=[num2str(RSPT_LAYER(k)) 'th layer'];
                end
                title(['Averaged from ' num2str(Pic_time_b(4)) ':00 ' mon(Pic_time_b(2),:) ' ' num2str(Pic_time_b(3)) ' ' num2str(Pic_time_b(1)) ' to ' num2str(Pic_time_e(4)) ':00 ' mon(Pic_time_e(2),:) ' ' num2str(Pic_time_e(3)) ' ' num2str(Pic_time_e(1)) ' (GMT +8)   ' RSPT_LAYER_NAME])
                hold on  
                
                %---------绘制岸线、岛------------------------
                fill(land(:,1),land(:,2),bgcolor);
                plot(land(:,1),land(:,2),'color',[0.2 0.2 0.2]);
                
                hold on
                
                for k2 =1:nisland;
                    eval(['fill( island' int2str(k2) '(:,1),island' int2str(k2) '(:,2),[' num2str(bgcolor) ']);']);
                    eval(['plot( island' int2str(k2) '(:,1),island' int2str(k2) '(:,2),''color'',[0.2 0.2 0.2]);']);
                    hold on
                end         
                
                fill(nanhuibiantan_weiken(:,1),nanhuibiantan_weiken(:,2),bgcolor);
                plot(nanhuibiantan_weiken(:,1),nanhuibiantan_weiken(:,2),'color',[0.2 0.2 0.2]);
                hold on
                
                if (FPT_COORDINATE=='BL') %按照图像的坐标轴选取绘制相应的label
                    xlabel('Longitude (\circE)');
                    ylabel('Latitude (\circN)');
                elseif (FPT_COORDINATE=='XY')
                    xlabel('Distance (km)');
                    ylabel('Distance (km)');
                end
                
                eval(['print(gcf,''-dpng'',''' fullfile(RSPT_Pic_Path,RSPT_Comp_Filename(1:end-4)) '_' num2str(RSPT_LAYER(k)) ''');']); %保存图像
                
            end
        end
    end 
    %-------- End Residual Salinity Field Distribution Drawing -------------    
    
    
    %------------- Sediment Field Distribution Drawing -------------------
    if (LOG_SED=='T') %如果需要画含沙量的平面图
        disp('Sediment Field Distribution Drawing ...')
        SED_Pic_Path=fullfile(OUT_DIRE,'field_distri\sediment'); %含沙量平面图存放路径
        [s,mess,messid]=mkdir(SED_Pic_Path); %建立含沙量平面图的存放文件夹
        EPT_Comp_Filepath=fullfile(IN_DIRE,'field_distri\elevation'); %模式的水位场计算结果文件路径
        SED_Comp_Filepath=fullfile(IN_DIRE,'field_distri\sediment'); %模式的含沙量场计算结果文件路径
        
        
        SED_TIME=str2num(SED_TIME); %将读入的含沙量绘图时刻的字符数组转化为浮点数组
        Num_SED_TIME=length(SED_TIME); %总共需要画几个时刻的图像
        NSED_LAYER=str2num(SED_LAYER); %将读入的含沙量绘图层数的字符数组转化为浮点数组
        Num_SED_LAYER=length(NSED_LAYER); %总共需要画几个层次的图像
        SED_contourscale=[-100 0 0.1 0.5 1 1.5 2 3 5]; %需要画的等值线
            
        
        
        for i=1:Num_SED_TIME %按照需要画图的数量进行循环
            
            disp(['  The ' num2str(SED_TIME(i)) ' Hours Drawing'])
            
            %------ Elevation Field Distribution Data Reading ,For Check Tideflat Grid ----------
            EPT_Comp_Filename_SN=SED_TIME(i)/(N_FPT/3600); %按照模式输出的时间间隔，计算需要画图时刻对应的序号
            EPT_Comp_Filename=sprintf('el_field_%06.6d',EPT_Comp_Filename_SN); %模式计算输出的水位文件
            SED_Comp_Filename_SN=SED_TIME(i)/(N_FPT/3600); %按照模式输出的时间间隔，计算需要画图时刻对应的序号
            SED_Comp_Filename=sprintf('sed_field_%06.6d',SED_Comp_Filename_SN); %模式计算输出的水位文件           
            
            if ~(exist(fullfile(EPT_Comp_Filepath,EPT_Comp_Filename))&&exist(fullfile(SED_Comp_Filepath,SED_Comp_Filename))) %判断水位场文件和含沙量场文件是否都存在，不存在则输出提示信息
                disp(['    ' EPT_Comp_Filename ' or ' SED_Comp_Filename ' dose not exist, drawing failed, please check the setting!']) %数据文件不存在，输出提示信息
            else %数据文件存在，则绘图
                eval(['fid_EPT_Comp = fopen(''' fullfile(EPT_Comp_Filepath,EPT_Comp_Filename) ''',''r'',''b'');']); %打开文件，由于该二进制文件是用“BIG_ENDIAN”封装，所以打开是选择参数‘b’
                brecord=fread(fid_EPT_Comp,1,'integer*4'); %读记录大小信息
                thour=fread(fid_EPT_Comp,1,'real*4');
                EPT_Comp_Data=fread(fid_EPT_Comp,IM*JM,'real*4');
                erecord=fread(fid_EPT_Comp,1,'integer*4'); %读记录大小信息
                
                if (brecord~=erecord) %如果一条记录的前后所记录的大小不一致，则退出程序
                    disp(['Error in reading file ' fullfile(EPT_Comp_Filepath,EPT_Comp_Filename) ', please check!'])
                    break
                end
                EPT_Comp_Data=reshape(EPT_Comp_Data,IM,JM); %将水位数据转化为矩阵形式
                fclose(fid_EPT_Comp); %关闭计算输出的水位文件
                %---------------- End Elevation Field Distribution Data Reading ----------------------
                
                
                
                %--------------------- Sediment Field Distribution Data Reading -----------------------
                eval(['fid_SED_Comp = fopen(''' fullfile(SED_Comp_Filepath,SED_Comp_Filename) ''',''r'',''b'');']); %打开文件，由于该二进制文件是用“BIG_ENDIAN”封装，所以打开是选择参数‘b’
                brecord=fread(fid_SED_Comp,1,'integer*4'); %读记录大小信息
                thour=fread(fid_SED_Comp,1,'real*4');
                SED_Comp_Data=fread(fid_SED_Comp,IM*JM*KB,'real*4');
                erecord=fread(fid_SED_Comp,1,'integer*4'); %读记录大小信息
                
                if (brecord~=erecord) %如果一条记录的前后所记录的大小不一致，则退出程序
                    disp(['Error in reading file ' fullfile(SED_Comp_Filepath,SED_Comp_Filename) ', please check!'])
                    break
                end
                SED_Comp_Data=reshape(SED_Comp_Data,IM,JM,KB); %将含沙量转化为矩阵形式
                fclose(fid_SED_Comp); %关闭计算输出的水位文件
                %------------------- End Sediment Field Distribution Data Reading -----------------------
                
                
                for i1=1:IM
                    for j1=1:JM
                        if (h(i1,j1)+EPT_Comp_Data(i1,j1)<DMIN||h(i1,j1)<-10) %将潮滩和陆地处的含沙量赋值为nan
                            SED_Comp_Data(i1,j1,:)=NaN;
                        end
                        
                        for k=1:KB-1
                            if SED_Comp_Data(i1,j1,k)<0
                                SED_Comp_Data(i1,j1,k)=0;
                            end
                        end
                    end
                end  
                
                
                
                Pic_time=datevec(Init_day+thour/24); %根据模式设置的起始时间和该数据记录的thour计算实际的时间，Pic_time包含 年、月、日、时、分、秒 信息
                
                
                for k=1:Num_SED_LAYER
                    close(figure(1));
                    disp(['    Layer ' num2str(SED_LAYER(k))])
                    if (FPT_COLOR_TYPE=='C') %根据绘图的颜色类型，选择颜色包文件
                        map=load ('Colormap\cm_sed2_C.dat');
                    elseif (FPT_COLOR_TYPE=='G')
                        map=load ('Colormap\cm_sed_G.dat');
                    end
                     map=[bgcolor;map]; %添加背景颜色，让含沙为nan的位置按背景色画
                    colormap(map);  %读入颜色包
                    
                    
                    set(gca,'box','on','FontName','times new roman','FontSize',12);
                    set(gca,'color',bgcolor);
                    set(gca,'xlim',[FPT_XMIN FPT_XMAX],'xtick',linspace(FPT_XMIN,FPT_XMAX,4));
                    set(gca,'ylim',[FPT_YMIN FPT_YMAX],'ytick',linspace(FPT_YMIN,FPT_YMAX,4));
                    set(gcf,'inverthardcopy','off'); %保存图片时按照设置的颜色而不自动调节
                    set(gcf,'color',[1 1 1]);
                    hold on                    
                    
                    warning off; %画contourf时会有很多警告信息，但不影响画图，所以关掉warning
                    eval(['contourf(xr,yr,SED_Comp_Data(:,:,' num2str(SED_LAYER(k)) '),SED_contourscale);']);
                    shading flat; 
                                
                    cmin=0;
                    cmax=4;
                     caxis([cmin-(cmax-cmin)/(length(map)-1) cmax]); %为了把背景色放到颜色包，而又不显示在colorbar里，‘/20’是由于颜色包里有20个颜色
                    if (FPT_COLOR_TYPE=='C') %如果画彩色图，则需要画colorbar
                        hc=colorbar;
                         set(hc,'ylim',[cmin cmax],'FontName','times new roman','FontSize',12); %设置colorbar的显示范围
                    end
  
      
                    
                    eval(['[cc,hh]=contour(xr,yr,SED_Comp_Data(:,:,' num2str(SED_LAYER(k)) '),SED_contourscale);']);
                    clabel(cc,hh,'FontSize',10,'FontName','Times new roman');
                    set(hh,'EdgeColor',[0.3 0.3 0.3],'linewidth',1) %设置等值线的颜色和粗细
                    
                    
                    
                    % 					if (SED_LAYER(k)==1)
                    % 						SED_LAYER_NAME='表';
                    % 					elseif (SED_LAYER(k)==5)
                    % 						SED_LAYER_NAME='底';
                    % 					else
                    % 						SED_LAYER_NAME=['第' num2str(SED_LAYER(k))];
                    % 					end
                    % 					title([num2str(Pic_time(1)) '年' num2str(Pic_time(2)) '月' num2str(Pic_time(3)) '日' num2str(Pic_time(4)) '时 ' SED_LAYER_NAME '层含沙量分布'])
                    if (SED_LAYER(k)==1)
                        SED_LAYER_NAME='Surface';
                    elseif (SED_LAYER(k)==KB-1)
                        SED_LAYER_NAME='Bottom';
                    elseif (SED_LAYER(k)==2)
                        SED_LAYER_NAME='2nd layer';
                    elseif (SED_LAYER(k)==3)
                        SED_LAYER_NAME='3rd layer';
                    else
                        SED_LAYER_NAME=[num2str(SED_LAYER(k)) 'th layer'];
                    end
                    title([num2str(Pic_time(4)) ':00 ' mon(Pic_time(2),:) ' ' num2str(Pic_time(3)) ' ' num2str(Pic_time(1)) ' (GMT +8)   ' SED_LAYER_NAME])
                    hold on  
                    
                    %---------绘制岸线、岛------------------------
                    fill(land(:,1),land(:,2),bgcolor);
                    plot(land(:,1),land(:,2),'color',[0.2 0.2 0.2]);
                    
                    hold on
                    
                    for k2 =1:nisland;
                        eval(['fill( island' int2str(k2) '(:,1),island' int2str(k2) '(:,2),[' num2str(bgcolor) ']);']);
                        eval(['plot( island' int2str(k2) '(:,1),island' int2str(k2) '(:,2),''color'',[0.2 0.2 0.2]);']);
                        hold on
                    end         
                    
                    fill(nanhuibiantan_weiken(:,1),nanhuibiantan_weiken(:,2),bgcolor);
                    plot(nanhuibiantan_weiken(:,1),nanhuibiantan_weiken(:,2),'color',[0.2 0.2 0.2]);
                    hold on
                    
                    if (FPT_COORDINATE=='BL') %按照图像的坐标轴选取绘制相应的label
                        xlabel('Longitude (\circE)');
                        ylabel('Latitude (\circN)');
                    elseif (FPT_COORDINATE=='XY')
                        xlabel('Distance (km)');
                        ylabel('Distance (km)');
                    end
                    
                    eval(['print(gcf,''-dpng'',''' fullfile(SED_Pic_Path,SED_Comp_Filename) '_' num2str(SED_LAYER(k)) ''');']); %保存图像
                    
                end
            end
        end      
    end 
    %------------- End Sediment Field Distribution Drawing ----------------   
    
    %----------- Residual Sedimnet Field Distribution Drawing -------------
      
    if (LOG_RSED=='T') %如果需要画含沙量的平面图
        
        disp('Residual sediment Field Distribution Drawing ...')
        RSED_Pic_Path=fullfile(OUT_DIRE,'residual_distri\field\sediment'); %余含沙量平面图存放路径
        [s,mess,messid]=mkdir(RSED_Pic_Path); %建立含沙量平面图的存放文件夹
        RSED_Comp_Filepath=fullfile(IN_DIRE,'resi_flux'); %模式的余含沙量场计算结果文件路径
        
 
        RSED_contourscale=[-0.1 0 0.1 0.5 1 1.5 2 5]; %需要画的等值线
        
        RSED_Comp_Fileinfo=dir(fullfile(RSED_Comp_Filepath,'resi_sed_flux_3d_*.out')); %模式计算结果文件信息（包括文件名，修改时间，大小，是否为文件路径）
        
        
        for i=1:length(RSED_Comp_Fileinfo) %按照需要画图的数量进行循环
            
            %--------------------- Residual sediment Field Distribution Data Reading -----------------------
            RSED_Comp_Filename=RSED_Comp_Fileinfo(i).name; %模式计算输出的余含沙量文件           
            eval(['fid_RSED_Comp = fopen(''' fullfile(RSED_Comp_Filepath,RSED_Comp_Filename) ''',''r'');']); %打开余含沙量文件
            RSED_Comp_Data=nan*ones(IM,JM,KB-1);
            Note=fscanf(fid_RSED_Comp,'%s',7);
            RSED_B=fscanf(fid_RSED_Comp,'%f',1);
            Note=fscanf(fid_RSED_Comp,'%s',2);
            RSED_E=fscanf(fid_RSED_Comp,'%f',1);
            Note=fgetl(fid_RSED_Comp);
            Note=fgetl(fid_RSED_Comp);
            disp(['  Residual Sediment Distribution from ' num2str(RSED_B) ' to ' num2str(RSED_E) ' Hour Drawing'])
            
            while(~feof(fid_RSED_Comp)) %判读是否读到文件尾部
                I=fscanf(fid_RSED_Comp,'%d',1); %读I
                J=fscanf(fid_RSED_Comp,'%d',1); %读J
                Note=fscanf(fid_RSED_Comp,'%f',2); %读X,Y
                RSED_Comp_Data(I,J,1:KB-1)=fscanf(fid_RSED_Comp,'%f',KB-1); %读X,Y位置kbm1层的含沙量值
                Note=fgetl(fid_RSED_Comp); %读掉末尾的回车
            end
            
            fclose(fid_RSED_Comp); %关闭计算输出的余含沙量文件
            %------------------- End Residual sediment Field Distribution Data Reading -----------------------
            
            for i1=1:IM
                for j1=1:JM
                    if (h(i1,j1)<RS_MIN_DEP) %小于RS_MIN_DEP位置不画
                        RSED_Comp_Data(i1,j1,:)=NaN; 
                        
                        for k=1:KB-1
                            if RSED_Comp_Data(i1,j1,k)<0
                                RSED_Comp_Data(i1,j1,k)=0;
                            end
                        end
                    end
                end
            end  
            
            
            
            
            Pic_time_b=datevec(Init_day+RSED_B/24); %根据模式设置的起始时间计算余含沙量场统计的起始时间
            Pic_time_e=datevec(Init_day+RSED_E/24); %根据模式设置的起始时间计算余含沙量场统计的结束时间
           
        RSED_LAYER=str2num(RSED_LAYER); %将读入的含沙量绘图层数的字符数组转化为浮点数组
        Num_RSED_LAYER=length(RSED_LAYER); %总共需要画几个层次的图像
            
            for k=1:Num_RSED_LAYER
                close(figure(1));
                disp(['    Layer ' num2str(RSED_LAYER(k))])
                if (FPT_COLOR_TYPE=='C') %根据绘图的颜色类型，选择颜色包文件
                    map=load ('Colormap\cm_sed2_C.dat');
                elseif (FPT_COLOR_TYPE=='G')
                    map=load ('Colormap\cm_sed_G.dat');
                end
                map=[bgcolor;map]; %添加背景颜色，让含沙量为nan的位置按背景色画
                colormap(map);  %读入颜色包
                
                set(gcf,'position',[50 50 1000 600]);
                x0=0.05;  y0=0.1; width=0.85; height=0.8;
                    
                 subplot('position',[x0 y0 width height]);
                 axis equal;
                set(gca,'box','on','FontName','times new roman','FontSize',12);
                set(gca,'color',bgcolor);
                set(gca,'xlim',[FPT_XMIN FPT_XMAX],'xtick',linspace(FPT_XMIN,FPT_XMAX,4));
                set(gca,'ylim',[FPT_YMIN FPT_YMAX],'ytick',linspace(FPT_YMIN,FPT_YMAX,4));
                set(gcf,'inverthardcopy','off'); %保存图片时按照设置的颜色而不自动调节
                set(gcf,'color',[1 1 1]);

                hold on 
                
                
                warning off; %画contourf时会有很多警告信息，但不影响画图，所以关掉warning
                eval(['contourf(xr,yr,RSED_Comp_Data(:,:,' num2str(RSED_LAYER(k)) '),RSED_contourscale);']);
                shading flat; 
%                 caxis([-0.1 5]); %为了把背景色放到颜色包，而又不显示在colorbar里，‘/20’是由于颜色包里有20个颜色
%                 if (FPT_COLOR_TYPE=='C') %如果画彩色图，则需要画colorbar
%                     hc=colorbar;
%                     set(hc,'ylim',[0 5],'FontName','times new roman','FontSize',12); %设置colorbar的显示范围
%                 end
                cmin=0;
                cmax=5;
                caxis([cmin-(cmax-cmin)/(length(map)-1) cmax]); %为了把背景色放到颜色包，而又不显示在colorbar里，‘/20’是由于颜色包里有20个颜色
                if (FPT_COLOR_TYPE=='C') %如果画彩色图，则需要画colorbar
                    hc=colorbar;
                    set(hc,'ylim',[cmin cmax],'FontName','times new roman','FontSize',12); %设置colorbar的显示范围
                end
                
                
                eval(['[cc,hh]=contour(xr,yr,RSED_Comp_Data(:,:,' num2str(RSED_LAYER(k)) '),RSED_contourscale);']);
                clabel(cc,hh,'FontSize',10,'FontName','Times new roman');
                set(hh,'EdgeColor',[0.3 0.3 0.3],'linewidth',1) %设置等值线的颜色和粗细
                
                
                
                if (RSED_LAYER(k)==1)
                    RSED_LAYER_NAME='Surface';
                elseif (RSED_LAYER(k)==KB-1)
                    RSED_LAYER_NAME='Bottom';
                elseif (RSED_LAYER(k)==2)
                    RSED_LAYER_NAME='2nd layer';
                elseif (RSED_LAYER(k)==3)
                    RSED_LAYER_NAME='3rd layer';
                else
                    RSED_LAYER_NAME=[num2str(RSED_LAYER(k)) 'th layer'];
                end
                title(['Averaged from ' num2str(Pic_time_b(4)) ':00 ' mon(Pic_time_b(2),:) ' ' num2str(Pic_time_b(3)) ' ' num2str(Pic_time_b(1)) ' to ' num2str(Pic_time_e(4)) ':00 ' mon(Pic_time_e(2),:) ' ' num2str(Pic_time_e(3)) ' ' num2str(Pic_time_e(1)) ' (GMT +8)   ' RSED_LAYER_NAME])
                hold on  
                
                %---------绘制岸线、岛------------------------
                fill(land(:,1),land(:,2),bgcolor);
                plot(land(:,1),land(:,2),'color',[0.2 0.2 0.2]);
                
                hold on
                
                for k2 =1:nisland;
                    eval(['fill( island' int2str(k2) '(:,1),island' int2str(k2) '(:,2),[' num2str(bgcolor) ']);']);
                    eval(['plot( island' int2str(k2) '(:,1),island' int2str(k2) '(:,2),''color'',[0.2 0.2 0.2]);']);
                    hold on
                end         
                
                fill(nanhuibiantan_weiken(:,1),nanhuibiantan_weiken(:,2),bgcolor);
                plot(nanhuibiantan_weiken(:,1),nanhuibiantan_weiken(:,2),'color',[0.2 0.2 0.2]);
                
                plot(deep_water_way(:,1),deep_water_way(:,2),'color',[0.2 0.2 0.2]);
                hold on
                
                if (FPT_COORDINATE=='BL') %按照图像的坐标轴选取绘制相应的label
                    xlabel('Longitude (\circE)');
                    ylabel('Latitude (\circN)');
                elseif (FPT_COORDINATE=='XY')
                    xlabel('Distance (km)');
                    ylabel('Distance (km)');
                end
                
                eval(['print(gcf,''-dpng'',''' fullfile(RSED_Pic_Path,RSED_Comp_Filename(1:end-4)) '_' num2str(RSED_LAYER(k)) ''');']); %保存图像
                
            end
        end
    end 
    %-------- End Residual sediment Field Distribution Drawing -------------       
    
    %------------- Bottom Shear Stress Field Distribution Drawing -------------------
    if (LOG_TAU=='T') %如果需要画切应力的平面图
        disp('Bottom Shear Stress Field Distribution Drawing ...')
        TAU_Pic_Path=fullfile(OUT_DIRE,'field_distri\stress'); %切应力平面图存放路径
        [s,mess,messid]=mkdir(TAU_Pic_Path); %建立切应力平面图的存放文件夹
        EPT_Comp_Filepath=fullfile(IN_DIRE,'field_distri\elevation'); %模式的水位场计算结果文件路径
        TAU_Comp_Filepath=fullfile(IN_DIRE,'field_distri\stress'); %模式的切应力场计算结果文件路径
        
        
        TAU_TIME=str2num(TAU_TIME); %将读入的切应力绘图时刻的字符数组转化为浮点数组
        Num_TAU_TIME=length(TAU_TIME); %总共需要画几个时刻的图像
        %        TAU_LAYER=str2num(TAU_LAYER); %将读入的切应力绘图层数的字符数组转化为浮点数组
        %        Num_TAU_LAYER=length(TAU_LAYER); %总共需要画几个层次的图像
        TAU_contourscale=[-0.1 0 0.5 1 3 5]; %需要画的等值线
        
        for i=1:Num_TAU_TIME %按照需要画图的数量进行循环
            
            disp(['  The ' num2str(TAU_TIME(i)) ' Hours Drawing'])
            
            %------ Elevation Field Distribution Data Reading ,For Check Tideflat Grid ----------
            EPT_Comp_Filename_SN=TAU_TIME(i)/(N_FPT/3600); %按照模式输出的时间间隔，计算需要画图时刻对应的序号
            EPT_Comp_Filename=sprintf('el_field_%06.6d',EPT_Comp_Filename_SN); %模式计算输出的水位文件
            TAU_Comp_Filename_SN=TAU_TIME(i)/(N_FPT/3600); %按照模式输出的时间间隔，计算需要画图时刻对应的序号
            TAU_Comp_Filename=sprintf('tau_field_%06.6d',TAU_Comp_Filename_SN); %模式计算输出的水位文件           
            
            if ~(exist(fullfile(EPT_Comp_Filepath,EPT_Comp_Filename))&&exist(fullfile(TAU_Comp_Filepath,TAU_Comp_Filename))) %判断水位场文件和盐度场文件是否都存在，不存在则输出提示信息
                disp(['    ' EPT_Comp_Filename ' or ' TAU_Comp_Filename ' dose not exist, drawing failed, please check the setting!']) %数据文件不存在，输出提示信息
            else %数据文件存在，则绘图
                eval(['fid_EPT_Comp = fopen(''' fullfile(EPT_Comp_Filepath,EPT_Comp_Filename) ''',''r'',''b'');']); %打开文件，由于该二进制文件是用“BIG_ENDIAN”封装，所以打开是选择参数‘b’
                brecord=fread(fid_EPT_Comp,1,'integer*4'); %读记录大小信息
                thour=fread(fid_EPT_Comp,1,'real*4');
                EPT_Comp_Data=fread(fid_EPT_Comp,IM*JM,'real*4');
                erecord=fread(fid_EPT_Comp,1,'integer*4'); %读记录大小信息
                
                if (brecord~=erecord) %如果一条记录的前后所记录的大小不一致，则退出程序
                    disp(['Error in reading file ' fullfile(EPT_Comp_Filepath,EPT_Comp_Filename) ', please check!'])
                    break
                end
                EPT_Comp_Data=reshape(EPT_Comp_Data,IM,JM); %将水位数据转化为矩阵形式
                fclose(fid_EPT_Comp); %关闭计算输出的水位文件
                %---------------- End Elevation Field Distribution Data Reading ----------------------
                
                
                
                %--------------------- TAUiment Field Distribution Data Reading -----------------------
                eval(['fid_TAU_Comp = fopen(''' fullfile(TAU_Comp_Filepath,TAU_Comp_Filename) ''',''r'',''b'');']); %打开文件，由于该二进制文件是用“BIG_ENDIAN”封装，所以打开是选择参数‘b’
                brecord=fread(fid_TAU_Comp,1,'integer*4'); %读记录大小信息
                thour=fread(fid_TAU_Comp,1,'real*4');
                TAU_Comp_Data=fread(fid_TAU_Comp,IM*JM,'real*4');
                erecord=fread(fid_TAU_Comp,1,'integer*4'); %读记录大小信息
                
                if (brecord~=erecord) %如果一条记录的前后所记录的大小不一致，则退出程序
                    disp(['Error in reading file ' fullfile(TAU_Comp_Filepath,TAU_Comp_Filename) ', please check!'])
                    break
                end
                TAU_Comp_Data=reshape(TAU_Comp_Data,IM,JM); %将切应力转化为矩阵形式
                fclose(fid_TAU_Comp); %关闭计算输出的水位文件
                %------------------- End Salinity Field Distribution Data Reading -----------------------
                
                
                for i1=1:IM
                    for j1=1:JM
                        if (h(i1,j1)+EPT_Comp_Data(i1,j1)<DMIN||h(i1,j1)<-10) %将潮滩和陆地处的切应力赋值为nan
                            TAU_Comp_Data(i1,j1,:)=NaN;
                        end
                    end
                end  
                
                
                
                Pic_time=datevec(Init_day+thour/24); %根据模式设置的起始时间和该数据记录的thour计算实际的时间，Pic_time包含 年、月、日、时、分、秒 信息
                
                
                %				for k=1:Num_TAU_LAYER
                close(figure(1));
                %					disp(['    Layer ' num2str(TAU_LAYER(k))])
                if (FPT_COLOR_TYPE=='C') %根据绘图的颜色类型，选择颜色包文件
                    map=load ('Colormap\cm_sed_C.dat');
                elseif (FPT_COLOR_TYPE=='G')
                    map=load ('Colormap\cm_sed_G.dat');
                end
                map=[bgcolor;map]; %添加背景颜色，让盐度为nan的位置按背景色画
                colormap(map);  %读入颜色包
                
                
                set(gca,'box','on','FontName','times new roman','FontSize',12);
                set(gca,'color',bgcolor);
                set(gca,'xlim',[FPT_XMIN FPT_XMAX],'xtick',linspace(FPT_XMIN,FPT_XMAX,4));
                set(gca,'ylim',[FPT_YMIN FPT_YMAX],'ytick',linspace(FPT_YMIN,FPT_YMAX,4));
                set(gcf,'inverthardcopy','off'); %保存图片时按照设置的颜色而不自动调节
                set(gcf,'color',[1 1 1]);
                hold on 
                
                
                warning off; %画contourf时会有很多警告信息，但不影响画图，所以关掉warning
                eval(['contourf(xr,yr,TAU_Comp_Data(:,:),TAU_contourscale);']);
                shading flat; 
                caxis([-0.1 5]); %为了把背景色放到颜色包，而又不显示在colorbar里，‘/20’是由于颜色包里有20个颜色
                if (FPT_COLOR_TYPE=='C') %如果画彩色图，则需要画colorbar
                    hc=colorbar;
                    set(hc,'ylim',[0 5],'FontName','times new roman','FontSize',12,'YColor',[0 0 0]); %设置colorbar的显示范围
                end
                
                
                eval(['[cc,hh]=contour(xr,yr,TAU_Comp_Data(:,:),TAU_contourscale);']);
                clabel(cc,hh,'FontSize',10,'FontName','Times new roman');
                set(hh,'EdgeColor',[0.3 0.3 0.3],'linewidth',1) %设置等值线的颜色和粗细
                

                hold on  
                
                %---------绘制岸线、岛------------------------
                fill(land(:,1),land(:,2),bgcolor);
                plot(land(:,1),land(:,2),'color',[0.2 0.2 0.2]);
                
                hold on
                
                for k2 =1:nisland;
                    eval(['fill( island' int2str(k2) '(:,1),island' int2str(k2) '(:,2),[' num2str(bgcolor) ']);']);
                    eval(['plot( island' int2str(k2) '(:,1),island' int2str(k2) '(:,2),''color'',[0.2 0.2 0.2]);']);
                    hold on
                end         
                
                fill(nanhuibiantan_weiken(:,1),nanhuibiantan_weiken(:,2),bgcolor);
                plot(nanhuibiantan_weiken(:,1),nanhuibiantan_weiken(:,2),'color',[0.2 0.2 0.2]);
                hold on
                
                if (FPT_COORDINATE=='BL') %按照图像的坐标轴选取绘制相应的label
                    xlabel('Longitude (\circE)');
                    ylabel('Latitude (\circN)');
                elseif (FPT_COORDINATE=='XY')
                    xlabel('Distance (km)');
                    ylabel('Distance (km)');
                end
                
                eval(['print(gcf,''-dpng'',''' fullfile(TAU_Pic_Path,TAU_Comp_Filename) ''');']); %保存图像
                
            end
        end      
    end
    %------------- End Bottom Shear Stress Field Distribution Drawing ----------------



