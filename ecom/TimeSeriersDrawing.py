%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%================== Time Seriers Drawing ====================

if (LOG_TSR_EL=='T'||LOG_TSR_VEL=='T'||LOG_TSR_S=='T') %站位时间序列绘图开关（T画，F不画）
    
    disp('================== Time Seriers Drawing ====================')
    
    TSR_Pic_Path=fullfile(OUT_DIRE,'timeseries'); %绘图结果存放路径
    
    [s,mess,messid]=mkdir(TSR_Pic_Path); %建立绘图结果存放文件夹
    
    TSR_Comp_Filepath=fullfile(IN_DIRE,'timeseries'); %模式计算结果文件路径
    
    TSR_Comp_Fileinfo=dir(fullfile(TSR_Comp_Filepath,'*.out')); %模式计算结果文件信息（包括文件名，修改时间，大小，是否为文件路径）
    
    TSR_Num=length(TSR_Comp_Fileinfo); %计算有几个待画图的文件
    
    TSR_LAYER_VEL=str2num(TSR_LAYER_VEL); %将读入的字符数组转化为浮点数数组
    
    TSR_LAYER_VEL=reshape(TSR_LAYER_VEL,1,length(TSR_LAYER_VEL)); %确保其为行向量
    
    TSR_LAYER_S=str2num(TSR_LAYER_S); %将读入的字符数组转化为浮点数数组
    
    TSR_LAYER_S=reshape(TSR_LAYER_S,1,length(TSR_LAYER_S)); %确保其为行向量
    
    for i=1:TSR_Num %按照文件个数进行循环画图
        
        TSR_Comp_Filename=TSR_Comp_Fileinfo(i).name; %计算结果文件的文件名
        
        Sitename=TSR_Comp_Filename(9:length(TSR_Comp_Filename)-4); %从文件名里取出站点名称
        
        disp(['Drawing Site ' Sitename ' ...'])
        
        eval(['fid_TSR_Comp = fopen(''' fullfile(TSR_Comp_Filepath,TSR_Comp_Filename) ''');']); %打开计算结果文件
        
        %   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        %         if(Sitename(2)=='n') %modified for the case
        %             TSR_LAG=53;
        %         elseif(Sitename(2)=='s')
        %             TSR_LAG=47;
        %         end
        %   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
        for i=1:4 %读文件的注释信息
            NOTE=fgetl(fid_TSR_Comp);
        end
        
        TSR_Compdata = fscanf(fid_TSR_Comp,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f',[17 inf]); %读出数据
        
        fclose(fid_TSR_Comp); %关闭模式计算结果文件
        
        
        %---------------------------Elevation Drawing----------------------------------    
        if (LOG_TSR_EL=='T')  %判断是否画水位过程线
            close(figure(1)); %关掉打开的绘图窗口
            disp('  Elevation Timeseries Drawing ...')
            subplot(312),
            plot((TSR_Compdata(1,:)/24)-TSR_LAG,TSR_Compdata(2,:),'k') %绘制计算的水位过程线图
            set(gca,'xlim',[TSR_BEG TSR_END],'xtick',linspace(TSR_BEG,TSR_END,TSR_END-TSR_BEG+1),'FontName','times new roman','FontSize',12)
            xlabel('Time (Day)');
            ylabel('Elevation (m)');
            %             eval(['title(''' Sitename '水位过程线'');'])
            eval(['title(''' Sitename ' Elevation Time Series'');'])
            
            hold on
            
            if (LOG_OBS_EL=='T') %水位验证点绘图开关（T画，F不画）
                
                TSR_Obs_Filepath='./Observation_Data/Elevation_Sequence'; %站位验证点存放路径，路径手动给出
                TSR_Obs_Filename=[Sitename '.dat']; %站位验证点文件名
                
                if exist(fullfile(TSR_Obs_Filepath,TSR_Obs_Filename)) %如果存在匹配的站位验证点数据则画出验证点，没有不画
                    eval(['fid_TSR_Obs = fopen(''' fullfile(TSR_Obs_Filepath,TSR_Obs_Filename) ''');']); %打开验证文件
                    TSR_Obsdata = fscanf(fid_TSR_Obs,'%f %f',[2 inf]); %读取验证点数据
                    fclose(fid_TSR_Obs); %关闭文件
                    plot(TSR_Obsdata(1,:)/24,TSR_Obsdata(2,:),'k.') %画出验证点
                    eval(['title(''' Sitename ' Elevation Validation'');'])
                end
                
            end
            
            eval(['print(gcf,''-dpng'',''' fullfile(TSR_Pic_Path,['El_' Sitename]) ''');'])
            
        end
        
        %---------------------------End Elevation Drawing-------------------------------    
             
        %---------------------------Velocity Drawing----------------------------------    
        if (LOG_TSR_VEL=='T') %判断是否画流速流向过程线
            close(figure(1));
            disp('  Velocity Timeseries Drawing ...')
            disp(['    LAYER   ' num2str(TSR_LAYER_VEL)])
            Num_TSR_LAYER_VEL=length(TSR_LAYER_VEL); %计算需要画几层
            posf = [350 36 600 900]; 
            set(gcf,'position',[posf(1) posf(2) posf(3) posf(4)]); %设置画布的位置和大小
            set(gcf,'paperpositionmode','auto');  %按照画布大小保存图形
            set(gcf,'inverthardcopy','on'); %保存图形时去掉绘图窗口的灰色
            
            
            for i=1:Num_TSR_LAYER_VEL %按照需要画图的层数进行循环
                eval(['subplot(' int2str(Num_TSR_LAYER_VEL*2) '1' int2str((i-1)*2+1) '),']);
                plot((TSR_Compdata(1,:)/24)-TSR_LAG,TSR_Compdata(2+(TSR_LAYER_VEL(i)-1)*3+1,:),'k') %绘制计算的第i层流速过程线图
                set(gca,'xlim',[TSR_BEG TSR_END],'xtick',linspace(TSR_BEG,TSR_END,TSR_END-TSR_BEG+1),'xticklabel',[],'FontName','times new roman','FontSize',12)
                ylabel('Speed (m/s)');
                
                if(i==1) %第一张图画出title
                    eval(['title(''' Sitename ' Velocity Time Series'');'])    
                end
                
                hold on
                
                
                eval(['subplot(' int2str(Num_TSR_LAYER_VEL*2) '1' int2str((i-1)*2+2) '),']);
                plot((TSR_Compdata(1,:)/24)-TSR_LAG,TSR_Compdata(2+(TSR_LAYER_VEL(i)-1)*3+2,:),'k') %绘制计算的第i层流向过程线图
                set(gca,'FontName','times new roman','FontSize',12)
                ylabel('Direction (\circ)');
                
                if(i==Num_TSR_LAYER_VEL )  %最后一张图画出x轴信息
                    set(gca,'xlim',[TSR_BEG TSR_END],'xtick',linspace(TSR_BEG,TSR_END,TSR_END-TSR_BEG+1),'ylim',[0 360],'ytick',linspace(0,360,5),'FontName','times new roman','FontSize',12)
                    xlabel('Time (Day)');        
                else
                    set(gca,'xlim',[TSR_BEG TSR_END],'xtick',linspace(TSR_BEG,TSR_END,TSR_END-TSR_BEG+1),'xticklabel',[],'ylim',[0 360],'ytick',linspace(0,360,5),'FontName','times new roman','FontSize',12)
                end
                
                hold on              
                
            end            
            
            
            if (LOG_OBS_VEL=='T') %流速流向验证点绘图开关（T画，F不画）
                
                TSR_Obs_Filepath='./Observation_Data/Velocity_Sequence'; %站位验证点存放路径，路径手动给出
                TSR_Obs_Filename=[Sitename '.dat']; %站位验证点文件名
                
                if exist(fullfile(TSR_Obs_Filepath,TSR_Obs_Filename)) %如果存在匹配的站位验证点数据则画出验证点，没有不画
                    %                     eval(['fid_TSR_Obs = fopen(''' fullfile(TSR_Obs_Filepath,TSR_Obs_Filename) ''');']); %打开验证文件
                    %                     TSR_Obsdata = fscanf(fid_TSR_Obs,'%f %f %f %f %f %f %f %f %f %f %f %f %f',[13 inf]); %读取验证点数据
                    %                     fclose(fid_TSR_Obs); %关闭文件
                    
                    TSR_Obsdata = load(fullfile(TSR_Obs_Filepath,TSR_Obs_Filename));
                    
                    for i=1:Num_TSR_LAYER_VEL %按照需要画图的层数进行循环
                        
                        eval(['subplot(' int2str(Num_TSR_LAYER_VEL*2) '1' int2str((i-1)*2+1) '),']);
                        %                         plot((TSR_Obsdata(1,:)/24),TSR_Obsdata(1+(TSR_LAYER_VEL(i)-1)*2+1,:),'b') %实测第i层流速数据
                        %                         plot((TSR_Obsdata(1,:)/24),TSR_Obsdata(1+(TSR_LAYER_VEL(i)-1)*2+3,:),'r') %实测第i+1层流速数据
                        plot((TSR_Obsdata(:,1)/24),TSR_Obsdata(:,1+(TSR_LAYER_VEL(i)-1)*2+1),'b') %实测第i层流速数据
                        plot((TSR_Obsdata(:,1)/24),TSR_Obsdata(:,1+(TSR_LAYER_VEL(i)-1)*2+3),'r') %实测第i+1层流速数据
                        
                        if(i==1) %第一张图画出title
                            eval(['title(''' Sitename ' Velocity Validation'');'])    
                        end
                        
                        hold on
                        
                        eval(['subplot(' int2str(Num_TSR_LAYER_VEL*2) '1' int2str((i-1)*2+2) '),']);
                        %                         plot((TSR_Obsdata(1,:)/24),TSR_Obsdata(1+(TSR_LAYER_VEL(i)-1)*2+2,:),'b') %实测第i层流向数据
                        %                         plot((TSR_Obsdata(1,:)/24),TSR_Obsdata(1+(TSR_LAYER_VEL(i)-1)*2+4,:),'r') %实测第i+1层流向数据
                        plot((TSR_Obsdata(:,1)/24),TSR_Obsdata(:,1+(TSR_LAYER_VEL(i)-1)*2+2),'b') %实测第i层流向数据
                        plot((TSR_Obsdata(:,1)/24),TSR_Obsdata(:,1+(TSR_LAYER_VEL(i)-1)*2+4),'r') %实测第i+1层流向数据
                        
                        
                        hold on              
                    end            
                    
                end
                
            end
            
            eval(['print(gcf,''-dpng'',''' fullfile(TSR_Pic_Path,['Vel_' Sitename]) ''');'])
            
        end
        %---------------------------End Veloci Drawing-------------------------------           
        
        
        
        
        
        
        %---------------------------Salinty Drawing----------------------------------    
        if (LOG_TSR_S=='T') %判断是否画盐度过程线
            close(figure(1))
            disp('  Salinity Timeseries Drawing ...')
            disp(['    LAYER   ' num2str(TSR_LAYER_S)])
            Num_TSR_LAYER_S=length(TSR_LAYER_S); %计算需要画几层
            posf = [350 36 600 900];
            set(gcf,'position',[posf(1) posf(2) posf(3) posf(4)]); %设置画布的位置和大小
            set(gcf,'paperpositionmode','auto');  %按照画布大小保存图形
            set(gcf,'inverthardcopy','on'); %保存图形时去掉绘图窗口的灰色
            
            for i=1:Num_TSR_LAYER_S %按照需要画图的层数进行循环
                eval(['subplot(' int2str(Num_TSR_LAYER_S) '1' int2str(i) '),']); 
                plot((TSR_Compdata(1,:)/24)-TSR_LAG,TSR_Compdata(2+(TSR_LAYER_S(i)-1)*3+3,:),'k') %绘制计算的盐度第i层过程线图
                set(gca,'FontName','times new roman','FontSize',12)
                if(i==1) %第一张图画出title
                    eval(['title(''' Sitename ' Salinity Time Series'');'])    
                end
                
                
                if(i==Num_TSR_LAYER_S )  %最后一张图画出x轴信息                 
                    set(gca,'xlim',[TSR_BEG TSR_END],'xtick',linspace(TSR_BEG,TSR_END,TSR_END-TSR_BEG+1),'FontName','times new roman','FontSize',12)
                    xlabel('Time (Day)');        
                else
                    
                    set(gca,'xlim',[TSR_BEG TSR_END],'xtick',linspace(TSR_BEG,TSR_END,TSR_END-TSR_BEG+1),'xticklabel',[],'FontName','times new roman','FontSize',12)
                end
                ylabel('Salinty');
                
                hold on
                
            end
            
            if (LOG_OBS_S=='T') %站位验证点绘图开关（T画，F不画）
                
                TSR_Obs_Filepath='./Observation_Data/Salinity_Sequence'; %站位验证点存放路径，路径手动给出
                TSR_Obs_Filename=[Sitename '.dat']; %站位验证点文件名
                
                if exist(fullfile(TSR_Obs_Filepath,TSR_Obs_Filename)) %如果存在匹配的站位验证点数据则画出验证点，没有不画
                    %                     eval(['fid_TSR_Obs = fopen(''' fullfile(TSR_Obs_Filepath,TSR_Obs_Filename) ''');']); %打开验证文件
                    %                     TSR_Obsdata = fscanf(fid_TSR_Obs,'%f %f %f %f %f %f %f',[7 inf]); %读取验证点数据
                    %                     fclose(fid_TSR_Obs); %关闭文件
                    
                    TSR_Obsdata = load(fullfile(TSR_Obs_Filepath,TSR_Obs_Filename));
                    
                    for i=1:Num_TSR_LAYER_S %按照需要画图的层数进行循环
                        
                        eval(['subplot(' int2str(Num_TSR_LAYER_S) '1' int2str((i)) '),']);
                        %                        plot((TSR_Obsdata(1,:)/24),TSR_Obsdata(1+TSR_LAYER_S(i),:),'b') %实测第i层盐度数据
                        %                        plot((TSR_Obsdata(1,:)/24),TSR_Obsdata(1+TSR_LAYER_S(i)+1,:),'r') %实测第i+1层盐度数据
                        plot((TSR_Obsdata(:,1)/24),TSR_Obsdata(:,1+TSR_LAYER_S(i)),'b') %实测第i层盐度数据
                        plot((TSR_Obsdata(:,1)/24),TSR_Obsdata(:,1+TSR_LAYER_S(i)+1),'r') %实测第i+1层盐度数据
                        set(gca,'FontName','times new roman','FontSize',12)
                        
                        if(i==1) %第一张图画出title
                            eval(['title(''' Sitename ' Salinity Validation'');'])    
                        end
                        
                        hold on
                        
                    end            
                    
                end
                
            end
            
            eval(['print(gcf,''-dpng'',''' fullfile(TSR_Pic_Path,['S_' Sitename]) ''');'])
            
        end
        %---------------------------End Salinity Drawing-------------------------------           
        
    end
    
    disp('================== End Of Time Seriers Drawing ====================')
    
end
%================== End Of Time Seriers Drawing ==================== 

if (LOG_TSR_SEC=='T') %断面时间序列绘图开关（T画，F不画）
    
    disp('================== Section Timeseries Drawing ====================')
    TSR_SEC_Pic_Path=fullfile(OUT_DIRE,'secflux'); %绘图结果存放路径
    [s,mess,messid]=mkdir(TSR_SEC_Pic_Path); %建立绘图结果存放文件夹
    TSR_SEC_Comp_Filepath=fullfile(IN_DIRE,'secflux'); %模式计算结果文件路径
    TSR_SEC_Comp_Fileinfo=dir(fullfile(TSR_SEC_Comp_Filepath,'*.out')); %模式计算结果文件信息（包括文件名，修改时间，大小，是否为文件路径）
    TSR_SEC_Num=length(TSR_SEC_Comp_Fileinfo); %计算有几个待画图的文件
    
    for i=1:TSR_SEC_Num %按照文件个数进行循环画图
        
        close(figure(1))
        TSR_SEC_Comp_Filename=TSR_SEC_Comp_Fileinfo(i).name; %计算结果文件的文件名
        Secname=TSR_SEC_Comp_Filename(9:length(TSR_SEC_Comp_Filename)-4); %从文件名里取出站点名称
        disp(['  Drawing Section ' Secname ' ...'])
        eval(['fid_TSR_SEC_Comp = fopen(''' fullfile(TSR_SEC_Comp_Filepath,TSR_SEC_Comp_Filename) ''');']); %打开计算结果文件
        
        for i=1:4 %读文件的注释信息
            NOTE=fgetl(fid_TSR_SEC_Comp);
        end
        
        TSR_SEC_Compdata = fscanf(fid_TSR_SEC_Comp,'%f %f %f %f %f',[5 inf]); %读出数据
        fclose(fid_TSR_SEC_Comp); %关闭模式计算结果文件
        
        posf = [350 36 600 900];
        set(gcf,'position',[posf(1) posf(2) posf(3) posf(4)]); %设置画布的位置和大小
        set(gcf,'paperpositionmode','auto');  %按照画布大小保存图形
        set(gcf,'inverthardcopy','on'); %保存图形时去掉绘图窗口的灰色
        hold on
        
        subplot(411)
        plot(TSR_SEC_Compdata(1,:)/24-TSR_SEC_LAG,TSR_SEC_Compdata(2,:),'k')
        set(gca,'xlim',[TSR_SEC_BEG TSR_SEC_END],'xtick',linspace(TSR_SEC_BEG,TSR_SEC_END,TSR_SEC_END-TSR_SEC_BEG+1),'xticklabel',[],'FontName','times new roman','FontSize',12)      
        ylabel('Flux (m^3/s)')
        title(['Section ' Secname ' Time Series'])
        
        subplot(412)
        plot(TSR_SEC_Compdata(1,:)/24-TSR_SEC_LAG,TSR_SEC_Compdata(3,:),'k')
        set(gca,'xlim',[TSR_SEC_BEG TSR_SEC_END],'xtick',linspace(TSR_SEC_BEG,TSR_SEC_END,TSR_SEC_END-TSR_SEC_BEG+1),'xticklabel',[],'FontName','times new roman','FontSize',12)      
        ylabel('Accumulated Flux (m^3)')
        
        subplot(413)
        plot(TSR_SEC_Compdata(1,:)/24-TSR_SEC_LAG,TSR_SEC_Compdata(4,:),'k')
        set(gca,'xlim',[TSR_SEC_BEG TSR_SEC_END],'xtick',linspace(TSR_SEC_BEG,TSR_SEC_END,TSR_SEC_END-TSR_SEC_BEG+1),'xticklabel',[],'FontName','times new roman','FontSize',12)       
        ylabel('Salt Flux (kg/s)')
        
        subplot(414)
        plot(TSR_SEC_Compdata(1,:)/24-TSR_SEC_LAG,TSR_SEC_Compdata(5,:),'k')
        set(gca,'xlim',[TSR_SEC_BEG TSR_SEC_END],'xtick',linspace(TSR_SEC_BEG,TSR_SEC_END,TSR_SEC_END-TSR_SEC_BEG+1),'FontName','times new roman','FontSize',12)      
        xlabel('Time (Day)');   
        ylabel('Accumulated Saltflux (kg)')
        
        eval(['print(gcf,''-dpng'',''' fullfile(TSR_SEC_Pic_Path,Secname) ''');'])
    end    
    disp('=============== End Of Section Time Seriers Drawing =================')
end


end

