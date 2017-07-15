% function [ output_args ] = SedSecDistribtionDrawing( input_args )
% 
% global IN_DIRE OUT_DIRE IYEAR IMONTH IDAY0 IM JM KB;
% global LOG_TSR_EL LOG_OBS_EL LOG_TSR_VEL LOG_OBS_VEL TSR_LAYER_VEL;
% global LOG_TSR_S LOG_OBS_S TSR_LAYER_S LOG_TSR_SEC TSR_BEG TSR_END TSR_LAG;
% global TSR_SEC_BEG TSR_SEC_END TSR_SEC_LAG N_FPT DMIN;
% global LOG_OPT LOG_EPT; 
% global EPT_TIME LOG_VPT_UV LOG_VPT_SD VPT_TIME VPT_LAYER VPT_INTERVAL;
% global VPT_SCALE 
% global LOG_FIELD LOG_SPT SPT_TIME SPT_LAYER LOG_RSPT RSPT_LAYER LOG_SED;
% global SED_TIME SED_LAYER LOG_RSED RSED_LAYER RS_MIN_DEP LOG_TAU TAU_TIME;
% global FPT_COORDINATE FPT_XMIN FPT_XMAX FPT_YMIN FPT_YMAX FPT_COLOR_TYPE;
% global LOG_SEC_S SEC_S_TIME LOG_SEC_RS LOG_SEC_SED SEC_SED_TIME LOG_SEC_RSED_FLUX;
% global SEC_P_TYPE SEC_NUM SEC_RESOLUTION SEC_MAX_HIGHT SEC_CONTROL_POINTS;
% global xr yr ; 

 if (LOG_SEC_RSED_FLUX=='T') %如果需要画余含沙量断面通量图
    %----------读取网格中心坐标和水深数据-------------------------
    Init_day=datenum([IYEAR IMONTH IDAY0]); %计算模式起算的天数       
    fid_XYH_BIN=fopen('ch_hzbc_griddepth'); %打开网格中心坐标、水深文件，该文件由模式输出
    brecord=fread(fid_XYH_BIN,1,'integer*4'); %读记录大小信息
    h=fread(fid_XYH_BIN,IM*JM,'real*4'); %读水深数据
    xr=fread(fid_XYH_BIN,IM*JM,'real*4'); %读网格中心点的X坐标
    yr=fread(fid_XYH_BIN,IM*JM,'real*4'); %读网格中心点的Y坐标
    
%    lon=fread(fid_XYH_BIN,IM*JM,'real*4'); %读网格中心点的X坐标
%    lat=fread(fid_XYH_BIN,IM*JM,'real*4'); %读网格中心点的Y坐标
    
    erecord=fread(fid_XYH_BIN,1,'integer*4'); %读记录大小信息
    fclose(fid_XYH_BIN);
    
    if (brecord~=erecord) %如果一条记录的前后所记录的大小不一致，则退出程序
        disp('Error in reading file ch_hzbc_griddepth, please check!')
%       break
    end
    

   if (FPT_COLOR_TYPE=='C') %按照绘图的颜色类型设置底图的背景色
        bgcolor=[0 1 0];
    elseif (FPT_COLOR_TYPE=='G')
        bgcolor=[0.3 0.3 0.3];
    end
     
    h=reshape(h,IM,JM); %转换成相应的矩阵形式
    xr=reshape(xr,IM,JM); %转换成相应的矩阵形式
    yr=reshape(yr,IM,JM); %转换成相应的矩阵形式
    
 

        xr=xr/1000; %将54坐标单位从m转换为km
        yr=yr/1000; %将54坐标单位从m转换为km
    
%----------------------END------------------------

        disp('Residual Sediment Flux Section Distribution Drawing ...')
        RSED_Flux_Sec_Pic_Path=fullfile(OUT_DIRE,'residual_distri\section\sediment_flux'); %余含沙量断面通量图存放路径（指定）
        [s,mess,messid]=mkdir(RSED_Flux_Sec_Pic_Path); %建立余含沙量断面图的存放文件夹（指定）
        
        REPT_Comp_Filepath=fullfile(IN_DIRE,'resi_euler'); %模式的余水位场计算结果文件路径
        REPT_Comp_Fileinfo=dir(fullfile(REPT_Comp_Filepath,'resi_euler_current_*.out')); %模式计算的余水位场（包括文件名，修改时间，大小，是否为文件路径）
        RSED_Flux_Comp_Filepath=fullfile(IN_DIRE,'resi_flux'); %模式的余含沙量断面通量计算结果文件路径
        RSED_Flux_Comp_Fileinfo=dir(fullfile(RSED_Flux_Comp_Filepath,'resi_sed_flux_3d_*.out')); %模式计算的含沙量断面通量（包括文件名，修改时间，大小，是否为文件路径）        
        RSED_Flux_contourscale=[-0.1 0 0.1 0.5 1 1.5 2 5]; %需要画的等值线        
        %计算网格边长
        
        for i=2:IM
            for j=2:JM
                h1(i,j)=sqrt((xr(i-1,j)-xr(i,j))^2+(yr(i-1,j)-yr(i,j))^2);
                h2(i,j)=sqrt((xr(i,j-1)-xr(i,j))^2+(yr(i,j-1)-yr(i,j))^2);
            end
        end
       
        for i=1:length(RSED_Flux_Comp_Fileinfo) %按照需要画图的数量进行循环
            
            %-----------------Residual Elevation Field Distribution Data Reading ------------------
            REPT_Comp_Filename=REPT_Comp_Fileinfo(i).name; %模式计算输出的余水位文件  
            eval(['fid_REPT_Comp = fopen(''' fullfile(REPT_Comp_Filepath,REPT_Comp_Filename) ''',''r'');']); %打开余水位文件
            REPT_Comp_Data=nan*ones(IM,JM);
            Note=fscanf(fid_REPT_Comp,'%s',7);
            REPT_B=fscanf(fid_REPT_Comp,'%f',1);
            Note=fscanf(fid_REPT_Comp,'%s',2);
            REPT_E=fscanf(fid_REPT_Comp,'%f',1);
            Note=fgetl(fid_REPT_Comp);
            Note=fgetl(fid_REPT_Comp);         
            
            while(~feof(fid_REPT_Comp)) %判读是否读到文件尾部
                I=fscanf(fid_REPT_Comp,'%d',1); %读I
                J=fscanf(fid_REPT_Comp,'%d',1); %读J
                Note=fscanf(fid_REPT_Comp,'%f',2); %读X,Y
                REPT_Comp_Data(I,J)=fscanf(fid_REPT_Comp,'%f',1); %读X,Y位置的水位值
                Note=fgetl(fid_REPT_Comp); %读掉其余部分
            end
            
            fclose(fid_REPT_Comp); %关闭计算输出的余水位文件
            %------------------- End Residal Elevation Field Distribution Data Reading -----------------------
           
           
            %-----------------Residual Sediment Field Distribution Data Reading ------------------
            RSED_Flux_Comp_Filename=RSED_Flux_Comp_Fileinfo(i).name; %模式计算输出的余含沙量文件  
            eval(['fid_RSED_Flux_Comp = fopen(''' fullfile(RSED_Flux_Comp_Filepath,RSED_Flux_Comp_Filename) ''',''r'');']); %打开余含沙量文件
            RSED_Flux_U_Comp_Data=nan*ones(IM,JM,KB-1);
            RSED_Flux_V_Comp_Data=nan*ones(IM,JM,KB-1);
            Note=fscanf(fid_RSED_Flux_Comp,'%s',7);
            RSED_Flux_B=fscanf(fid_RSED_Flux_Comp,'%f',1);
            Note=fscanf(fid_RSED_Flux_Comp,'%s',2);
            RSED_Flux_E=fscanf(fid_RSED_Flux_Comp,'%f',1);
            Note=fgetl(fid_RSED_Flux_Comp);
            Note=fgetl(fid_RSED_Flux_Comp);
            disp(['  Residual Sediment Section Distribution from ' num2str(RSED_Flux_B) ' to ' num2str(RSED_Flux_E) ' Hour Drawing'])
            
            while(~feof(fid_RSED_Flux_Comp)) %判读是否读到文件尾部
                I=fscanf(fid_RSED_Flux_Comp,'%d',1); %读I
                J=fscanf(fid_RSED_Flux_Comp,'%d',1); %读J
                Note=fscanf(fid_RSED_Flux_Comp,'%f',2+KB-1); %读X,Y和余含沙量数据
                for kk=1:KB-1
                    RSED_Flux_U_Comp_Data(I,J,kk)=fscanf(fid_RSED_Flux_Comp,'%f',1); %读X,Y位置U方向的余含沙通量
                    RSED_Flux_V_Comp_Data(I,J,1:kk)=fscanf(fid_RSED_Flux_Comp,'%f',1); %读X,Y位置V方向的余含沙通量
                end
                Note=fgetl(fid_RSED_Flux_Comp); %读掉末尾的回车
            end
            
            fclose(fid_RSED_Flux_Comp); %关闭计算输出的余含沙量文件

            %------------------- End Residual Sediment Field Distribution Data Reading -----------------------
            
            fsm=ones(IM,JM);
            for i1=1:IM
                for j1=1:JM
                    if isnan(RSED_Flux_U_Comp_Data(i1,j1))==1 %按照fsmadd判断
                        fsm(i1,j1)=0; 
                        REPT_Comp_Data(i1,j1)=-SEC_MAX_HIGHT; 
                        RSED_Flux_U_Comp_Data(i1,j1,:)=NaN;    
                        RSED_Flux_V_Comp_Data(i1,j1,:)=NaN;    
                    end
                    
                    if (h(i1,j1)<RS_MIN_DEP) %小于RS_MIN_DEP位置不画
                        fsm(i1,j1)=0;
                        REPT_Comp_Data(i1,j1)=-SEC_MAX_HIGHT;
                        RSED_Flux_U_Comp_Data(i1,j1,:)=NaN;    
                        RSED_Flux_V_Comp_Data(i1,j1,:)=NaN;    
                    end
                end
            end  
        
            Pic_time_b=datevec(Init_day+RSED_Flux_B/24); %根据模式设置的起始时间计算余含沙量场统计的起始时间
            Pic_time_e=datevec(Init_day+RSED_Flux_E/24); %根据模式设置的起始时间计算余含沙量场统计的结束时间
            

            close(figure(2));
  
            
            %=================== 绘制断面图 ======================
             
                 for nsec=1:SEC_NUM %按照断面数目进行循环
                    expression=['SEC_CONTROL_POINT=str2num(SEC_CONTROL_POINTS{',num2str(nsec),',1});'];%将断面控制点统一赋值到SEC_CONTROL_POINT中
                    eval(expression); 
                    
%                   expression=['SEC_CONTROL_POINT=SEC_CONTROL_POINTS' num2str(nsec) ';'] %将断面控制点统一赋值到SEC_CONTROL_POINT中
%                   eval(expression);                    
%                    figure(1),plot(SEC_CONTROL_POINT(:,1),SEC_CONTROL_POINT(:,2),'linewidth',1.5,'color','k'); %在图形1中（含沙量表层平面图）中画出断面的位置
%                    hold on
%                    text(SEC_CONTROL_POINT(1,1),SEC_CONTROL_POINT(1,2),num2str(nsec),'FontName','times new roman','FontSize',15);
%                    hold on
                    
                    SEC_SED_POSSITION_X=ones(0); %初始化断面（根据分辨率将控制点扩充）横坐标向量，初值为空
                    SEC_SED_POSSITION_Y=ones(0); %初始化断面（根据分辨率将控制点扩充）纵坐标向量，初值为空
                    for ii=1:length(SEC_CONTROL_POINT)-1 %对断面相邻两个控制点进行处理
                        xs=SEC_CONTROL_POINT(ii,1); %起点横坐标
                        xe=SEC_CONTROL_POINT(ii+1,1); %终点横坐标
                        ys=SEC_CONTROL_POINT(ii,2); %起点纵坐标
                        ye=SEC_CONTROL_POINT(ii+1,2); %终点纵坐标
                        dis=sqrt((xs-xe)^2+(ys-ye)^2); %起点和终点的距离
                        num=ceil(dis/SEC_RESOLUTION); %按照预先设定的分辨率确定每两个控制点间分几段
                        dx=(xe-xs)/num; %x轴方向每段的长度
                        dy=(ye-ys)/num; %y轴方向每段的长度
                        
                        for jj=1:num
                            xi=xs+(jj-1)*dx; %计算扩充点的横坐标
                            yi=ys+(jj-1)*dy; %计算扩充点的纵坐标
                            SEC_SED_POSSITION_X=[SEC_SED_POSSITION_X xi]; %添加到横坐标向量
                            SEC_SED_POSSITION_Y=[SEC_SED_POSSITION_Y yi]; %添加到纵坐标向量
                        end
                        
                    end
                    SEC_SED_POSSITION_X=[SEC_SED_POSSITION_X SEC_CONTROL_POINT(end,1)]; %将最后一个控制点添加到横坐标向量
                    SEC_SED_POSSITION_Y=[SEC_SED_POSSITION_Y SEC_CONTROL_POINT(end,2)]; %将最后一个控制点添加到纵坐标向量

                    sec_h=plain_interp(SEC_SED_POSSITION_X,SEC_SED_POSSITION_Y,h,IM,JM,KB,xr,yr,fsm,h,h1,h2,SEC_MAX_HIGHT); %对水深进行双线性插值，陆地处用预先设定的输出图像最高高程定义（SEC_MAX_HIGHT）
                    sec_el=plain_interp(SEC_SED_POSSITION_X,SEC_SED_POSSITION_Y,REPT_Comp_Data,IM,JM,KB,xr,yr,fsm,h,h1,h2,-SEC_MAX_HIGHT); %对水位进行双线性插值，陆地处用预先设定的输出图像最高高程的相反数定义（-SEC_MAX_HIGHT） 注：考虑到水位与水深的正方向相反
                    sec_sed_flux_u=plain_interp(SEC_SED_POSSITION_X,SEC_SED_POSSITION_Y,RSED_Flux_U_Comp_Data,IM,JM,KB,xr,yr,fsm,h,h1,h2,nan); %对余含沙通量进行双线性插值，陆地处用nan表示
                    sec_sed_flux_v=plain_interp(SEC_SED_POSSITION_X,SEC_SED_POSSITION_Y,RSED_Flux_V_Comp_Data,IM,JM,KB,xr,yr,fsm,h,h1,h2,nan); %对余含沙通量进行双线性插值，陆地处用nan表示
                    
                    %注：由于含沙量定义在网格的中央，为了得到完整的图像，垂向上上下都等值外插半个网格
                    for kk=1:length(SEC_SED_POSSITION_X)
                        distance(kk)=sqrt((SEC_SED_POSSITION_X(kk)-SEC_SED_POSSITION_X(1))^2+(SEC_SED_POSSITION_Y(kk)-SEC_SED_POSSITION_Y(1))^2); %计算断面各点距离断面起点的距离
                    end
                    x_sec_s=ones(KB+1,1)*distance; %有效的含沙量共有KB-1个，加上向两端外插的2个数值，所以垂向为KB+1个数据
                 
                 
                    for ii=1:length(SEC_SED_POSSITION_X)
                       sec_d=sec_h(ii)+sec_el(ii); %总水深等于水深加上水位
                                       
                        for kk=1:KB-1    
                            y_sec_s(kk,ii)=-sec_el(ii)+sec_d/(KB-1)*(kk-1)+0.5*sec_d/(KB-1); %求中间KB-1个有效点的所在的位置（深度）
                        end
                    end

                    y_sec_s=[-sec_el;y_sec_s;sec_h]; %表层位置取水位的相反数，底层位置取相应的水深                    
                    sec_sed_flux_u=[sec_sed_flux_u(1,:);sec_sed_flux_u;sec_sed_flux_u(KB-1,:)]; %含沙量等值外插半个网格
                    sec_sed_flux_v=[sec_sed_flux_v(1,:);sec_sed_flux_v;sec_sed_flux_v(KB-1,:)]; %含沙量等值外插半个网格                    

                    %法向：站在断面起点，面向断面终点，右手为正
                    ddx=SEC_CONTROL_POINT(2,1)-SEC_CONTROL_POINT(1,1); %断面方向
                    ddy=SEC_CONTROL_POINT(2,2)-SEC_CONTROL_POINT(1,2); %断面方向
                    
                    for kk=1:KB+1
                        sec_sed_flux_t(kk,1)=(sec_sed_flux_u(kk,1)*ddx+sec_sed_flux_v(kk,1)*ddy)/sqrt(ddx^2+ddy^2); %断面切向沙通量
                        sec_sed_flux_n(kk,1)=(sec_sed_flux_v(kk,1)*ddx-sec_sed_flux_u(kk,1)*ddy)/sqrt(ddx^2+ddy^2); %断面法向沙通量
                    end
                    
                    for ii=2:length(SEC_SED_POSSITION_X)
                        ddx=SEC_SED_POSSITION_X(ii)-SEC_SED_POSSITION_X(ii-1); %断面方向
                        ddy=SEC_SED_POSSITION_Y(ii)-SEC_SED_POSSITION_Y(ii-1); %断面方向

                        for kk=1:KB+1
                            sec_sed_flux_t(kk,ii)=(sec_sed_flux_u(kk,ii)*ddx+sec_sed_flux_v(kk,ii)*ddy)/sqrt(ddx^2+ddy^2); %断面切向沙通量
                            sec_sed_flux_n(kk,ii)=(sec_sed_flux_u(kk,ii)*ddy-sec_sed_flux_v(kk,ii)*ddx)/sqrt(ddx^2+ddy^2); %断面法向沙通量
                        end
                    end
                    
                 
                    figure(2),clf;
                    sec_sed_map=load ('Colormap\cm_el_C.dat');
                    sec_sed_map=[[0 1 0];sec_sed_map];
                    colormap(sec_sed_map);
                   
                    set(gcf,'position',[50 50 1000 600]);
                    x0=0.08; y0=0.08; width=0.85;height=0.6;                    
                    subplot('position',[x0 y0 width height]);
    %               subplot(2,2,1:2);
                    
%                   contourf(x_sec_s,y_sec_s,sec_sed,contourscale); 
                    contourf(x_sec_s,y_sec_s,sec_sed_flux_n); %画法向流速等值线图
                    shading flat;
                    hold on;
                    Data_Max=max(max(max(sec_sed_flux_n))); %求出水位数据的最大值最小值
                    Data_Min=min(min(min(sec_sed_flux_n)));
                    cmax=ceil(max(abs(Data_Max),abs(Data_Min))); %求出最大的振幅，并向上取整
                    cmin=-cmax; %因为颜色包里的颜色是对称的，所以要使正负坐标对称
                    caxis([cmin-(cmax-cmin)/(length(sec_sed_map)-1) cmax]); %为了把背景色放到颜色包，而又不显示在colorbar里
                    hc=colorbar;
                    set(hc,'ylim',[cmin cmax],'FontName','times new roman','FontSize',12); %设置colorbar的显示范围
                    
                    set(gca,'color',[0 1 0],'FontName','times new roman','FontSize',12); %背景色设为与地形一样的绿色          
                    set(gca,'xlim',[0 max(max(x_sec_s))]); %横坐标设置为0至最大距离
                    set(gca,'ylim',[SEC_MAX_HIGHT ceil(max(max(y_sec_s)))]); %纵坐标设置为最高高程至水深最大值的向上取整
%                   set(gca,'ylim',[0 1]); %纵坐标设置为最高高程至水深最大值的向上取整
                    set(gca,'YDir','reverse') %将y轴反向
                    set(gcf,'inverthardcopy','off'); %保存图片时按照设置的颜色而不自动调节
                    set(gcf,'color',[1 1 1]); 
                    hold on;
                    
                 
                    sec_sed_flux_w=sec_sed_flux_t; 
                    for ii=1:length(SEC_SED_POSSITION_X)
                        for kk=1:KB+1
                            if isnan(sec_sed_flux_t(kk,ii))==0
                                sec_sed_flux_w(kk,ii)=0;%由于没有数据输出，暂时将垂向分量定义为0
                            end
                        end
                    end
                    
                    SEC_XMIN=0;
                    SEC_XMAX=max(max(x_sec_s));
                    SEC_YMIN=SEC_MAX_HIGHT;
                    SEC_YMAX=ceil(max(max(y_sec_s)));
                    
                    %在quiver函数中设定箭头按实际大小画，而不是自动调节大小，所以需要把将原始数据放大到和坐标轴相应的大小尺度
                    Draw_Flux_t=sec_sed_flux_t*(SEC_XMAX-SEC_XMIN)*0.03/VPT_SCALE; %(SEC_XMAX-SEC_XMIN)*0.03/VPT_SCALE表示将x轴的3%的长度对应VPT_SCALE的大小
                    Draw_Flux_w=sec_sed_flux_w*(SEC_YMAX-SEC_YMIN)*0.05/VPT_SCALE; %(SEC_YMAX-SEC_YMIN)*0.2/VPT_SCALE表示将y轴的5%的长度对应VPT_SCALE的大小    
                    plotvector(x_sec_s,y_sec_s,Draw_Flux_t,Draw_Flux_w,VPT_SCALE,1,'r','line','uv','fix',1); 
                    hold on
                    
                    fill([distance(1) distance distance(end)],[SEC_MAX_HIGHT -sec_el SEC_MAX_HIGHT],'w') %将水面以上填充为白色
                    hold on                 
                    %                     plot(distance,-sec_el,'r')
                    hold on
                    
                    fill([distance(1) distance distance(end)],[ceil(max(max(y_sec_s))) sec_h ceil(max(max(y_sec_s)))],'g') %将地形包括陆地填充为绿色
                    hold on
                    plot(distance,sec_h,'g') %填充时，轮廓线与填充色不一致，故把地形轮廓画为绿色用于遮盖
                    hold on

                    
                    x_scale=SEC_XMIN+(SEC_XMAX-SEC_XMIN)*0.03; %横向为从左向右3%的位置
                    y_scale=SEC_YMIN+(SEC_YMAX-SEC_YMIN)*0.92; %纵向为从上而下80%的位置
%                    quiver([x_scale x_scale],[y_scale y_scale],[(SEC_XMAX-SEC_XMIN)*0.03 0],[0 -(SEC_YMAX-SEC_YMIN)*0.05],0,'k'); %按照画箭头尺度（这里画图片纵横各5%）
                     quiver(x_scale,y_scale,(SEC_XMAX-SEC_XMIN)*0.03,0,0,'k'); %按照画箭头尺度（这里画图片纵横各3%） 只画水平方向的尺标
%                    plotvector(x_scale,y_scale,(SEC_XMAX-SEC_XMIN)*0.03,0,1,'k'); 
                    hold on
                    text(SEC_XMIN+(SEC_XMAX-SEC_XMIN)*0.03,SEC_YMIN+(SEC_YMAX-SEC_YMIN)*0.95,[num2str(1./VPT_SCALE) 'kg/m^2s'],'FontName','times new roman','FontSize',8); %标注箭头的量值
                    hold on
                    
                    xlabel('Distance (km)');
                    ylabel('Depth (m)')
%                       title([num2str(Pic_time(1)) '年' num2str(Pic_time(2)) '月' num2str(Pic_time(3)) '日' num2str(Pic_time(4)) '时 断面' num2str(nsec) '含沙量分布'])
%                    title(['Residual sediment flux during time from ' num2str(Pic_time_b(4)) ':00 ' mon(Pic_time_b(2),:) ' ' num2str(Pic_time_b(3)) ' ' num2str(Pic_time_b(1)) ' to ' num2str(Pic_time_e(4)) ':00 ' mon(Pic_time_e(2),:) ' ' num2str(Pic_time_e(3)) ' ' num2str(Pic_time_e(1)) ' Sec' num2str(nsec)])
                    
               figure(2),eval(['print(gcf,''-dpng'',''' fullfile(RSED_Flux_Sec_Pic_Path,RSED_Flux_Comp_Filename(1:end-4)) '_specify_sec' num2str(nsec) ''');']); %保存断面图像
               figure(2),eval(['saveas(gcf,''' fullfile(RSED_Flux_Sec_Pic_Path,RSED_Flux_Comp_Filename(1:end-4)) '_specify_sec' num2str(nsec) '.eps'',''psc2 '');']);  %保存断面图像
                    clear SEC_SED_POSSITION_X SEC_SED_POSSITION_Y x_sec_s y_sec_s distance sec_h sec_el sec_sed sec_sed_flux_t sec_sed_flux_n %清除部分变量
                  end
         end      
 end 
            close(figure(2));
            