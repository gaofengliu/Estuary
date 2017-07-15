
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
%------------- Salinity Section Distribution Drawing -------------------
    if (LOG_SEC_S=='T') %如果需要画盐度的断面图
        disp('Salinity Section Distribution Drawing ...')
        S_Sec_Pic_Path=fullfile(OUT_DIRE,'section_distri\salinity'); %盐度断面图存放路径（指定）
        [s,mess,messid]=mkdir(S_Sec_Pic_Path); %建立盐度断面图的存放文件夹（指定）
        EPT_Comp_Filepath=fullfile(IN_DIRE,'field_distri\elevation'); %模式的水位场计算结果文件路径
        SPT_Comp_Filepath=fullfile(IN_DIRE,'field_distri\salinity'); %模式的盐度场计算结果文件路径
        
        SEC_S_TIME=str2num(SEC_S_TIME); %将读入的盐度断面绘图时刻的字符数组转化为浮点数组
        Num_SEC_S_TIME=length(SEC_S_TIME); %总共需要画几个时刻的图像
        SPT_contourscale=[-0.5 0 0.5 1 3 5 10 15 20 25 30 35]; %需要画的等值线
        
        %计算网格边长
        for i=2:IM
            for j=2:JM
                h1(i,j)=sqrt((xr(i-1,j)-xr(i,j))^2+(yr(i-1,j)-yr(i,j))^2);
                h2(i,j)=sqrt((xr(i,j-1)-xr(i,j))^2+(yr(i,j-1)-yr(i,j))^2);
            end
        end
       
        if strcmp(SEC_P_TYPE,'SPECIFY') %如果已在设置文件中指定断面的位置
            for nsec=1:SEC_NUM
                eval(['SEC_CONTROL_POINTS_' num2str(nsec) '=str2num(SEC_CONTROL_POINTS_' num2str(nsec) ');']); %将读入的断面位置从字符变量改为浮点变量
            end
        end
        
        for i=1:Num_SEC_S_TIME %按照需要画图的数量进行循环
            
            disp(['  The ' num2str(SEC_S_TIME(i)) ' Hours Drawing'])
            
            %------ Elevation Field Distribution Data Reading ,For Check Tideflat Grid ----------
            EPT_Comp_Filename_SN=SEC_S_TIME(i)/(N_FPT/3600); %按照模式输出的时间间隔，计算需要画图时刻对应的序号
            EPT_Comp_Filename=sprintf('el_field_%06.6d',EPT_Comp_Filename_SN); %模式计算输出的水位文件
            SPT_Comp_Filename_SN=SEC_S_TIME(i)/(N_FPT/3600); %按照模式输出的时间间隔，计算需要画图时刻对应的序号
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
                
                fsm=ones(IM,JM);
                for i1=1:IM
                    for j1=1:JM
                        if (h(i1,j1)+EPT_Comp_Data(i1,j1)<DMIN||h(i1,j1)<-10) %将潮滩和陆地处的盐度赋值为nan和干湿网格标记
                            SPT_Comp_Data(i1,j1,:)=NaN;
                            fsm(i1,j1)=0;    
                        end
                        
                        for k=1:KB-1
                            if SPT_Comp_Data(i1,j1,k)<0
                                SPT_Comp_Data(i1,j1,k)=0;
                            end
                        end
                        
                    end
                end  
                
                
                
                Pic_time=datevec(Init_day+thour/24); %根据模式设置的起始时间和该数据记录的thour计算实际的时间，Pic_time包含 年、月、日、时、分、秒 信息
                
                
                close(figure(1));
                close(figure(2));
                figure(1),
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
                contourf(xr,yr,SPT_Comp_Data(:,:,1),SPT_contourscale);
                shading flat; 
                caxis([-0.5 35]); %为了把背景色放到颜色包，而又不显示在colorbar里，‘/20’是由于颜色包里有20个颜色
                if (FPT_COLOR_TYPE=='C') %如果画彩色图，则需要画colorbar
                    hc=colorbar;
                    set(hc,'ylim',[0 35],'FontName','times new roman','FontSize',12); %设置colorbar的显示范围
                end
                
                [cc,hh]=contour(xr,yr,SPT_Comp_Data(:,:,1),SPT_contourscale);
                clabel(cc,hh,'FontSize',10,'FontName','Times new roman');
                set(hh,'EdgeColor',[0.3 0.3 0.3],'linewidth',1) %设置等值线的颜色和粗细
                
                %  				title([num2str(Pic_time(1)) '年' num2str(Pic_time(2)) '月' num2str(Pic_time(3)) '日' num2str(Pic_time(4)) '时 表层盐度分布'])
                title([num2str(Pic_time(4)) ':00 ' mon(Pic_time(2),:) ' ' num2str(Pic_time(3)) ' ' num2str(Pic_time(1)) ' (GMT +8)   Surface'])
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
             
                %=================== 绘制断面图 ======================
                if strcmp(SEC_P_TYPE,'SPECIFY') %如果断面位置预先指定 
                    
                    for nsec=1:SEC_NUM %按照断面数目进行循环
                        eval(['SEC_CONTROL_POINTS=SEC_CONTROL_POINTS_' num2str(nsec) ';']) %将断面控制点统一赋值到SEC_CONTROL_POINTS中
                        figure(1),plot(SEC_CONTROL_POINTS(:,1),SEC_CONTROL_POINTS(:,2),'linewidth',1.5,'color','k'); %在图形1中（盐度表层平面图）中画出断面的位置
                        hold on
                        text(SEC_CONTROL_POINTS(1,1),SEC_CONTROL_POINTS(1,2),num2str(nsec),'FontName','times new roman','FontSize',15);
                        hold on
                        
                        SEC_S_POSSITION_X=ones(0); %初始化断面（根据分辨率将控制点扩充）横坐标向量，初值为空
                        SEC_S_POSSITION_Y=ones(0); %初始化断面（根据分辨率将控制点扩充）纵坐标向量，初值为空
                        for ii=1:length(SEC_CONTROL_POINTS)-1 %对断面相邻两个控制点进行处理
                            xs=SEC_CONTROL_POINTS(ii,1); %起点横坐标
                            xe=SEC_CONTROL_POINTS(ii+1,1); %终点横坐标
                            ys=SEC_CONTROL_POINTS(ii,2); %起点纵坐标
                            ye=SEC_CONTROL_POINTS(ii+1,2); %终点纵坐标
                            dis=sqrt((xs-xe)^2+(ys-ye)^2); %起点和终点的距离
                            num=ceil(dis/SEC_RESOLUTION); %按照预先设定的分辨率确定每两个控制点间分几段
                            dx=(xe-xs)/num; %x轴方向每段的长度
                            dy=(ye-ys)/num; %y轴方向每段的长度
                            
                            for jj=1:num
                                xi=xs+(jj-1)*dx; %计算扩充点的横坐标
                                yi=ys+(jj-1)*dy; %计算扩充点的纵坐标
                                SEC_S_POSSITION_X=[SEC_S_POSSITION_X xi]; %添加到横坐标向量
                                SEC_S_POSSITION_Y=[SEC_S_POSSITION_Y yi]; %添加到纵坐标向量
                            end
                            
                        end
                        SEC_S_POSSITION_X=[SEC_S_POSSITION_X SEC_CONTROL_POINTS(end,1)]; %将最后一个控制点添加到横坐标向量
                        SEC_S_POSSITION_Y=[SEC_S_POSSITION_Y SEC_CONTROL_POINTS(end,2)]; %将最后一个控制点添加到纵坐标向量
                        
                        sec_h=plain_interp(SEC_S_POSSITION_X,SEC_S_POSSITION_Y,h,IM,JM,KB,xr,yr,fsm,h,h1,h2,SEC_MAX_HIGHT); %对水深进行双线性插值，陆地处用预先设定的输出图像最高高程定义（SEC_MAX_HIGHT）
                        sec_el=plain_interp(SEC_S_POSSITION_X,SEC_S_POSSITION_Y,EPT_Comp_Data,IM,JM,KB,xr,yr,fsm,h,h1,h2,-SEC_MAX_HIGHT); %对水位进行双线性插值，陆地处用预先设定的输出图像最高高程的相反数定义（-SEC_MAX_HIGHT） 注：考虑到水位与水深的正方向相反
                        sec_s=plain_interp(SEC_S_POSSITION_X,SEC_S_POSSITION_Y,SPT_Comp_Data,IM,JM,KB,xr,yr,fsm,h,h1,h2,nan); %对盐度行双线性插值，陆地处用nan表示
                        
                        
                        %注：由于盐度定义在网格的中央，为了得到完整的图像，垂向上上下都等值外插半个网格
                        for kk=1:length(SEC_S_POSSITION_X)
                            distance(kk)=sqrt((SEC_S_POSSITION_X(kk)-SEC_S_POSSITION_X(1))^2+(SEC_S_POSSITION_Y(kk)-SEC_S_POSSITION_Y(1))^2); %计算断面各点距离断面起点的距离
                        end
                        x_sec_s=ones(KB+1,1)*distance; %有效的盐度空有KB-1个，加上向两端外插的2个数值，所以垂向为7个数据
                        
                        for ii=1:length(SEC_S_POSSITION_X)
                            sec_d=sec_h(ii)+sec_el(ii); %总水深等于水深加上水位
                            for kk=1:KB-1    
                                y_sec_s(kk,ii)=-sec_el(ii)+sec_d/(KB-1)*(kk-1)+0.5*sec_d/(KB-1); %求中间KB-1个有效点的所在的位置（深度）
                            end
                        end
                        y_sec_s=[-sec_el;y_sec_s;sec_h]; %表层位置取水位的相反数，底层位置取相应的水深
                        
                        
                        sec_s=[sec_s(1,:);sec_s;sec_s(KB-1,:)]; %盐度等值外插半个网格
                        
                        
                        %此处定义所需要的等值线范围，应保证有0.5等值线（淡水），还未考虑周全，因此目前画等值线时按照默认！       
                        if max(max(sec_s))<=1
                            contourscale=linspace(0,1,6);
                        elseif min(min(sec_s))<=0.5&&max(max(sec_s))>1
                            contourscale=[0 linspace(0.5,ceil(max(max(sec_s))),5)];
                        else
                            contourscale=linspace(floor(min(min(sec_s))),ceil(max(max(sec_s))),6);
                        end
                        
                        figure(2),clf;
                        sec_s_map=load ('Colormap\cm_s_C.dat');
                        sec_s_map=[[0 1 0];sec_s_map];
                        colormap(sec_s_map);
                        
                        subplot(3,1,2);
                        
                        %                        contourf(x_sec_s,y_sec_s,sec_s,contourscale); 
                        contourf(x_sec_s,y_sec_s,sec_s); 
                        shading flat;
                        hold on;
                        caxis([-0.5,35]); 
                        hc=colorbar;
                        set(hc,'ylim',[0 35],'FontName','times new roman','FontSize',12);
                        
                        set(gca,'color',[0 1 0],'FontName','times new roman','FontSize',12); %背景色设为与地形一样的绿色          
                        set(gca,'xlim',[0 max(max(x_sec_s))]); %横坐标设置为0至最大距离
                        set(gca,'ylim',[SEC_MAX_HIGHT ceil(max(max(y_sec_s)))]); %纵坐标设置为最高高程至水深最大值的向上取整
                        set(gca,'YDir','reverse') %将y轴反向
                        set(gcf,'inverthardcopy','off'); %保存图片时按照设置的颜色而不自动调节
                        set(gcf,'color',[1 1 1]); 
                        hold on;
                        
                        %                          [cc,hh]=contour(x_sec_s,y_sec_s,sec_s,contourscale,'k');
                        
                        [cc,hh]=contour(x_sec_s,y_sec_s,sec_s,'k');
                        hold on;
                        clabel(cc,hh,'fontsize',10,'fontname','Times new roman','color','k','rotation',0);
                        hold on;
                        
                        fill([distance(1) distance distance(end)],[SEC_MAX_HIGHT -sec_el SEC_MAX_HIGHT],'w') %将水面以上填充为白色
                        hold on                 
                        %                         plot(distance,-sec_el,'w')
                        hold on
                        
                        fill([distance(1) distance distance(end)],[ceil(max(max(y_sec_s))) sec_h ceil(max(max(y_sec_s)))],'g') %将地形包括陆地填充为绿色
                        hold on
                        plot(distance,sec_h,'g') %填充时，轮廓线与填充色不一致，故把地形轮廓画为绿色用于遮盖
                        hold on
                        xlabel('Distance (km)');
                        ylabel('Depth (m)')
                        %                       title([num2str(Pic_time(1)) '年' num2str(Pic_time(2)) '月' num2str(Pic_time(3)) '日' num2str(Pic_time(4)) '时 断面' num2str(nsec) '盐度分布'])
                        title([num2str(Pic_time(4)) ':00 ' mon(Pic_time(2),:) ' ' num2str(Pic_time(3)) ' ' num2str(Pic_time(1)) ' (GMT +8)   Sec' num2str(nsec)])                         
                        
                        figure(2),eval(['print(gcf,''-dpng'',''' fullfile(S_Sec_Pic_Path,SPT_Comp_Filename) '_specify_sec' num2str(nsec) ''');']); %保存断面图像
                        
                        clear SEC_S_POSSITION_X SEC_S_POSSITION_Y x_sec_s y_sec_s distance sec_h sec_el sec_s %清除部分变量
                        
                    end
                    figure(1),eval(['print(gcf,''-dpng'',''' fullfile(S_Sec_Pic_Path,SPT_Comp_Filename) '_specify_sec_check'');']); %保存平面图像
                    
             end
        end      
    end 
    %------------- End Salinity Section Distribution Drawing ----------------  
