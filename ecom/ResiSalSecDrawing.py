function [ output_args ] = ResiSalSecDrawing( input_args )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
 if (LOG_SEC_RS=='T') %如果需要画盐度的断面图
        disp('Residual Salinity Section Distribution Drawing ...')
        RS_Sec_Pic_Path_Manual=fullfile(OUT_DIRE,'residual_distri\section\salinity\manual'); %余盐度断面图存放路径（手动）
        [s,mess,messid]=mkdir(RS_Sec_Pic_Path_Manual); %建立余盐度断面图的存放文件夹（手动）
        RS_Sec_Pic_Path_Specify=fullfile(OUT_DIRE,'residual_distri\section\salinity\specify'); %余盐度断面图存放路径（指定）
        [s,mess,messid]=mkdir(RS_Sec_Pic_Path_Specify); %建立余盐度断面图的存放文件夹（指定）
        
        REPT_Comp_Filepath=fullfile(IN_DIRE,'resi_euler'); %模式的余水位场计算结果文件路径
        REPT_Comp_Fileinfo=dir(fullfile(REPT_Comp_Filepath,'resi_euler_current_*.out')); %模式计算的余水位场（包括文件名，修改时间，大小，是否为文件路径）
        RSPT_Comp_Filepath=fullfile(IN_DIRE,'resi_flux'); %模式的余盐度场计算结果文件路径
        RSPT_Comp_Fileinfo=dir(fullfile(RSPT_Comp_Filepath,'resi_salt_flux_3d_*.out')); %模式计算的余盐度场（包括文件名，修改时间，大小，是否为文件路径）
        
        
        RSPT_contourscale=[-0.5 0 0.5 1 3 5 10 15 20 25 30 35]; %需要画的等值线
        
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
        
        for i=1:length(RSPT_Comp_Fileinfo) %按照需要画图的数量进行循环
            
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
            
            %-----------------Residual Salinity Field Distribution Data Reading ------------------
            RSPT_Comp_Filename=RSPT_Comp_Fileinfo(i).name; %模式计算输出的余盐度文件  
            eval(['fid_RSPT_Comp = fopen(''' fullfile(RSPT_Comp_Filepath,RSPT_Comp_Filename) ''',''r'');']); %打开余盐度文件
            RSPT_Comp_Data=nan*ones(IM,JM,KB-1);
            Note=fscanf(fid_RSPT_Comp,'%s',7);
            RSPT_B=fscanf(fid_RSPT_Comp,'%f',1);
            Note=fscanf(fid_RSPT_Comp,'%s',2);
            RSPT_E=fscanf(fid_RSPT_Comp,'%f',1);
            Note=fgetl(fid_RSPT_Comp);
            Note=fgetl(fid_RSPT_Comp);
            disp(['  Residual Salinty Section Distribution from ' num2str(RSPT_B) ' to ' num2str(RSPT_E) ' Hour Drawing'])
            
            while(~feof(fid_RSPT_Comp)) %判读是否读到文件尾部
                I=fscanf(fid_RSPT_Comp,'%d',1); %读I
                J=fscanf(fid_RSPT_Comp,'%d',1); %读J
                Note=fscanf(fid_RSPT_Comp,'%f',2); %读X,Y
                RSPT_Comp_Data(I,J,1:KB-1)=fscanf(fid_RSPT_Comp,'%f',KB-1); %读X,Y位置kbm1层的盐度值
                Note=fgetl(fid_RSPT_Comp); %读掉末尾的回车
            end
            
            fclose(fid_RSPT_Comp); %关闭计算输出的余盐度文件
            %------------------- End Residual Salinity Field Distribution Data Reading -----------------------
            
            fsm=ones(IM,JM);
            for i1=1:IM
                for j1=1:JM
                    if (isnan(RSPT_Comp_Data(i1,j1))==1) %按照fsmadd判断
                        fsm(i1,j1)=0; 
                        REPT_Comp_Data(i1,j1)=NaN; 
                        RSPT_Comp_Data(i1,j1,:)=NaN;    
                    end
                    
                    
                    for k=1:KB-1
                        if RSPT_Comp_Data(i1,j1,k)<0
                            RSPT_Comp_Data(i1,j1,k)=0;
                        end
                    end
                    
                    
                    if (h(i1,j1)<RS_MIN_DEP) %小于RS_MIN_DEP位置不画
                        fsm(i1,j1)=0;
                        REPT_Comp_Data(i1,j1)=NaN;
                        RSPT_Comp_Data(i1,j1,:)=NaN; 
                    end
                    
                end
            end  
            
            
            
            
            Pic_time_b=datevec(Init_day+RSPT_B/24); %根据模式设置的起始时间计算余盐度场统计的起始时间
            Pic_time_e=datevec(Init_day+RSPT_E/24); %根据模式设置的起始时间计算余盐度场统计的结束时间
            
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
            contourf(xr,yr,RSPT_Comp_Data(:,:,1),RSPT_contourscale);
            shading flat; 
            caxis([-0.5 35]); %为了把背景色放到颜色包，而又不显示在colorbar里，‘/20’是由于颜色包里有20个颜色
            if (FPT_COLOR_TYPE=='C') %如果画彩色图，则需要画colorbar
                hc=colorbar;
                set(hc,'ylim',[0 35],'FontName','times new roman','FontSize',12); %设置colorbar的显示范围
            end
            
            [cc,hh]=contour(xr,yr,RSPT_Comp_Data(:,:,1),RSPT_contourscale);
            clabel(cc,hh,'FontSize',10,'FontName','Times new roman');
            set(hh,'EdgeColor',[0.3 0.3 0.3],'linewidth',1) %设置等值线的颜色和粗细
            
            %  				title([num2str(Pic_time(1)) '年' num2str(Pic_time(2)) '月' num2str(Pic_time(3)) '日' num2str(Pic_time(4)) '时 表层盐度分布'])
            title(['Averaged from ' num2str(Pic_time_b(4)) ':00 ' mon(Pic_time_b(2),:) ' ' num2str(Pic_time_b(3)) ' ' num2str(Pic_time_b(1)) ' to ' num2str(Pic_time_e(4)) ':00 ' mon(Pic_time_e(2),:) ' ' num2str(Pic_time_e(3)) ' ' num2str(Pic_time_e(1)) ' (GMT +8)   Surface'])
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
                    sec_el=plain_interp(SEC_S_POSSITION_X,SEC_S_POSSITION_Y,REPT_Comp_Data,IM,JM,KB,xr,yr,fsm,h,h1,h2,-SEC_MAX_HIGHT); %对水位进行双线性插值，陆地处用预先设定的输出图像最高高程的相反数定义（-SEC_MAX_HIGHT） 注：考虑到水位与水深的正方向相反
                    sec_s=plain_interp(SEC_S_POSSITION_X,SEC_S_POSSITION_Y,RSPT_Comp_Data,IM,JM,KB,xr,yr,fsm,h,h1,h2,nan); %对盐度行双线性插值，陆地处用nan表示
                    
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
%                     clabel(cc,hh,'fontsize',10,'fontname','Times new roman','color','k','rotation',0);
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
                    title(['Averaged from ' num2str(Pic_time_b(4)) ':00 ' mon(Pic_time_b(2),:) ' ' num2str(Pic_time_b(3)) ' ' num2str(Pic_time_b(1)) ' to ' num2str(Pic_time_e(4)) ':00 ' mon(Pic_time_e(2),:) ' ' num2str(Pic_time_e(3)) ' ' num2str(Pic_time_e(1)) ' (GMT +8)   Sec' num2str(nsec)])
                    
                    figure(2),eval(['print(gcf,''-dpng'',''' fullfile(RS_Sec_Pic_Path_Specify,RSPT_Comp_Filename(1:end-4)) '_specify_sec' num2str(nsec) ''');']); %保存断面图像
                    
                    clear SEC_S_POSSITION_X SEC_S_POSSITION_Y x_sec_s y_sec_s distance sec_h sec_el sec_s %清除部分变量
                    
                end
                figure(1),eval(['print(gcf,''-dpng'',''' fullfile(RS_Sec_Pic_Path_Specify,RSPT_Comp_Filename(1:end-4)) '_specify_sec_check'');']); %保存平面图像
                
            elseif strcmp(SEC_P_TYPE,'MANUAL')  
                nsec=0;
                fid_SEC_CONTROL_POINTS_MANUAL=fopen(fullfile(RS_Sec_Pic_Path_Manual,'Sec_control_points.dat'),'w'); %打开断面控制点数据记录文件
                while (1)
                    figure(1),SEC_CONTROL_POINTS=ginput(100); %鼠标点击断面位置，点击100次自动生成，不足100次按Enter结束
                    if (isempty(SEC_CONTROL_POINTS)~=1&&min(size(SEC_CONTROL_POINTS))>1) %点击2次以上则进行绘图
                        nsec=nsec+1; %断面图片编号
                        
                        fprintf(fid_SEC_CONTROL_POINTS_MANUAL,'%10.2f %10.2f;',SEC_CONTROL_POINTS'); %按照格式输入数据
                        fprintf(fid_SEC_CONTROL_POINTS_MANUAL,'\n'); %输入回车
                        
                        figure(1),plot(SEC_CONTROL_POINTS(:,1),SEC_CONTROL_POINTS(:,2),'linewidth',1.5,'color','k'); %在图形1中（盐度表层平面图）中画出断面的位置
                        hold on
                        text(SEC_CONTROL_POINTS(1,1),SEC_CONTROL_POINTS(1,2),num2str(nsec),'FontSize',15);
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
                        sec_el=plain_interp(SEC_S_POSSITION_X,SEC_S_POSSITION_Y,REPT_Comp_Data,IM,JM,KB,xr,yr,fsm,h,h1,h2,-SEC_MAX_HIGHT); %对水位进行双线性插值，陆地处用预先设定的输出图像最高高程的相反数定义（-SEC_MAX_HIGHT） 注：考虑到水位与水深的正方向相反
                        sec_s=plain_interp(SEC_S_POSSITION_X,SEC_S_POSSITION_Y,RSPT_Comp_Data,IM,JM,KB,xr,yr,fsm,h,h1,h2,nan); %对盐度行双线性插值，陆地处用nan表示
                        
                        
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
                        
                        %                           contourf(x_sec_s,y_sec_s,sec_s,contourscale); 
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
                        
                        %                            [cc,hh]=contour(x_sec_s,y_sec_s,sec_s,contourscale,'k');
                        
                        [cc,hh]=contour(x_sec_s,y_sec_s,sec_s,'k');
                        hold on;
                        clabel(cc,hh,'FontSize',10,'FontName','Times new roman','color','k','rotation',0);
                        hold on;
                        
                        fill([distance(1) distance distance(end)],[SEC_MAX_HIGHT -sec_el SEC_MAX_HIGHT],'w') %将水面以上填充为白色
                        hold on                 
                        %                           plot(distance,-sec_el,'w')
                        hold on
                        
                        fill([distance(1) distance distance(end)],[ceil(max(max(y_sec_s))) sec_h ceil(max(max(y_sec_s)))],'g') %将地形包括陆地填充为绿色
                        hold on
                        plot(distance,sec_h,'g') %填充时，轮廓线与填充色不一致，故把地形轮廓画为绿色用于遮盖
                        hold on
                        xlabel('Distance (km)');
                        ylabel('Depth (m)')
                        %  							title([num2str(Pic_time(1)) '年' num2str(Pic_time(2)) '月' num2str(Pic_time(3)) '日' num2str(Pic_time(4)) '时 断面' num2str(nsec) '盐度分布'])
                        title(['Averaged from ' num2str(Pic_time_b(4)) ':00 ' mon(Pic_time_b(2),:) ' ' num2str(Pic_time_b(3)) ' ' num2str(Pic_time_b(1)) ' to ' num2str(Pic_time_e(4)) ':00 ' mon(Pic_time_e(2),:) ' ' num2str(Pic_time_e(3)) ' ' num2str(Pic_time_e(1)) ' (GMT +8)   Sec' num2str(nsec)])
                        
                        figure(2),eval(['print(gcf,''-dpng'',''' fullfile(RS_Sec_Pic_Path_Manual,RSPT_Comp_Filename(1:end-4)) '_manual_sec' num2str(nsec) ''');']); %保存断面图像
                        
                        clear SEC_S_POSSITION_X SEC_S_POSSITION_Y x_sec_s y_sec_s distance sec_h sec_el sec_s %清除部分变量                             
                        
                    elseif (isempty(SEC_CONTROL_POINTS)~=1&&min(size(SEC_CONTROL_POINTS))==1) %点击1次屏显提示信息，重试
                        disp('    Single point can not make section, two point at least, please try again!')
                        
                    else %直接按Enter键退出
                        break
                    end
                end
                fclose(fid_SEC_CONTROL_POINTS_MANUAL); %关闭断面控制点数据记录文件
                figure(1),eval(['print(gcf,''-dpng'',''' fullfile(RS_Sec_Pic_Path_Manual,RSPT_Comp_Filename(1:end-4)) '_manual_sec_check'');']); %保存平面图像
            end
            
        end      
    end 
    %------------- End Residual Salinity Section Distribution Drawing ---------------- 

end

