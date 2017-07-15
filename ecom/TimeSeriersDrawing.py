%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%================== Time Seriers Drawing ====================

if (LOG_TSR_EL=='T'||LOG_TSR_VEL=='T'||LOG_TSR_S=='T') %վλʱ�����л�ͼ���أ�T����F������
    
    disp('================== Time Seriers Drawing ====================')
    
    TSR_Pic_Path=fullfile(OUT_DIRE,'timeseries'); %��ͼ������·��
    
    [s,mess,messid]=mkdir(TSR_Pic_Path); %������ͼ�������ļ���
    
    TSR_Comp_Filepath=fullfile(IN_DIRE,'timeseries'); %ģʽ�������ļ�·��
    
    TSR_Comp_Fileinfo=dir(fullfile(TSR_Comp_Filepath,'*.out')); %ģʽ�������ļ���Ϣ�������ļ������޸�ʱ�䣬��С���Ƿ�Ϊ�ļ�·����
    
    TSR_Num=length(TSR_Comp_Fileinfo); %�����м�������ͼ���ļ�
    
    TSR_LAYER_VEL=str2num(TSR_LAYER_VEL); %��������ַ�����ת��Ϊ����������
    
    TSR_LAYER_VEL=reshape(TSR_LAYER_VEL,1,length(TSR_LAYER_VEL)); %ȷ����Ϊ������
    
    TSR_LAYER_S=str2num(TSR_LAYER_S); %��������ַ�����ת��Ϊ����������
    
    TSR_LAYER_S=reshape(TSR_LAYER_S,1,length(TSR_LAYER_S)); %ȷ����Ϊ������
    
    for i=1:TSR_Num %�����ļ���������ѭ����ͼ
        
        TSR_Comp_Filename=TSR_Comp_Fileinfo(i).name; %�������ļ����ļ���
        
        Sitename=TSR_Comp_Filename(9:length(TSR_Comp_Filename)-4); %���ļ�����ȡ��վ������
        
        disp(['Drawing Site ' Sitename ' ...'])
        
        eval(['fid_TSR_Comp = fopen(''' fullfile(TSR_Comp_Filepath,TSR_Comp_Filename) ''');']); %�򿪼������ļ�
        
        %   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        %         if(Sitename(2)=='n') %modified for the case
        %             TSR_LAG=53;
        %         elseif(Sitename(2)=='s')
        %             TSR_LAG=47;
        %         end
        %   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
        for i=1:4 %���ļ���ע����Ϣ
            NOTE=fgetl(fid_TSR_Comp);
        end
        
        TSR_Compdata = fscanf(fid_TSR_Comp,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f',[17 inf]); %��������
        
        fclose(fid_TSR_Comp); %�ر�ģʽ�������ļ�
        
        
        %---------------------------Elevation Drawing----------------------------------    
        if (LOG_TSR_EL=='T')  %�ж��Ƿ�ˮλ������
            close(figure(1)); %�ص��򿪵Ļ�ͼ����
            disp('  Elevation Timeseries Drawing ...')
            subplot(312),
            plot((TSR_Compdata(1,:)/24)-TSR_LAG,TSR_Compdata(2,:),'k') %���Ƽ����ˮλ������ͼ
            set(gca,'xlim',[TSR_BEG TSR_END],'xtick',linspace(TSR_BEG,TSR_END,TSR_END-TSR_BEG+1),'FontName','times new roman','FontSize',12)
            xlabel('Time (Day)');
            ylabel('Elevation (m)');
            %             eval(['title(''' Sitename 'ˮλ������'');'])
            eval(['title(''' Sitename ' Elevation Time Series'');'])
            
            hold on
            
            if (LOG_OBS_EL=='T') %ˮλ��֤���ͼ���أ�T����F������
                
                TSR_Obs_Filepath='./Observation_Data/Elevation_Sequence'; %վλ��֤����·����·���ֶ�����
                TSR_Obs_Filename=[Sitename '.dat']; %վλ��֤���ļ���
                
                if exist(fullfile(TSR_Obs_Filepath,TSR_Obs_Filename)) %�������ƥ���վλ��֤�������򻭳���֤�㣬û�в���
                    eval(['fid_TSR_Obs = fopen(''' fullfile(TSR_Obs_Filepath,TSR_Obs_Filename) ''');']); %����֤�ļ�
                    TSR_Obsdata = fscanf(fid_TSR_Obs,'%f %f',[2 inf]); %��ȡ��֤������
                    fclose(fid_TSR_Obs); %�ر��ļ�
                    plot(TSR_Obsdata(1,:)/24,TSR_Obsdata(2,:),'k.') %������֤��
                    eval(['title(''' Sitename ' Elevation Validation'');'])
                end
                
            end
            
            eval(['print(gcf,''-dpng'',''' fullfile(TSR_Pic_Path,['El_' Sitename]) ''');'])
            
        end
        
        %---------------------------End Elevation Drawing-------------------------------    
             
        %---------------------------Velocity Drawing----------------------------------    
        if (LOG_TSR_VEL=='T') %�ж��Ƿ��������������
            close(figure(1));
            disp('  Velocity Timeseries Drawing ...')
            disp(['    LAYER   ' num2str(TSR_LAYER_VEL)])
            Num_TSR_LAYER_VEL=length(TSR_LAYER_VEL); %������Ҫ������
            posf = [350 36 600 900]; 
            set(gcf,'position',[posf(1) posf(2) posf(3) posf(4)]); %���û�����λ�úʹ�С
            set(gcf,'paperpositionmode','auto');  %���ջ�����С����ͼ��
            set(gcf,'inverthardcopy','on'); %����ͼ��ʱȥ����ͼ���ڵĻ�ɫ
            
            
            for i=1:Num_TSR_LAYER_VEL %������Ҫ��ͼ�Ĳ�������ѭ��
                eval(['subplot(' int2str(Num_TSR_LAYER_VEL*2) '1' int2str((i-1)*2+1) '),']);
                plot((TSR_Compdata(1,:)/24)-TSR_LAG,TSR_Compdata(2+(TSR_LAYER_VEL(i)-1)*3+1,:),'k') %���Ƽ���ĵ�i�����ٹ�����ͼ
                set(gca,'xlim',[TSR_BEG TSR_END],'xtick',linspace(TSR_BEG,TSR_END,TSR_END-TSR_BEG+1),'xticklabel',[],'FontName','times new roman','FontSize',12)
                ylabel('Speed (m/s)');
                
                if(i==1) %��һ��ͼ����title
                    eval(['title(''' Sitename ' Velocity Time Series'');'])    
                end
                
                hold on
                
                
                eval(['subplot(' int2str(Num_TSR_LAYER_VEL*2) '1' int2str((i-1)*2+2) '),']);
                plot((TSR_Compdata(1,:)/24)-TSR_LAG,TSR_Compdata(2+(TSR_LAYER_VEL(i)-1)*3+2,:),'k') %���Ƽ���ĵ�i�����������ͼ
                set(gca,'FontName','times new roman','FontSize',12)
                ylabel('Direction (\circ)');
                
                if(i==Num_TSR_LAYER_VEL )  %���һ��ͼ����x����Ϣ
                    set(gca,'xlim',[TSR_BEG TSR_END],'xtick',linspace(TSR_BEG,TSR_END,TSR_END-TSR_BEG+1),'ylim',[0 360],'ytick',linspace(0,360,5),'FontName','times new roman','FontSize',12)
                    xlabel('Time (Day)');        
                else
                    set(gca,'xlim',[TSR_BEG TSR_END],'xtick',linspace(TSR_BEG,TSR_END,TSR_END-TSR_BEG+1),'xticklabel',[],'ylim',[0 360],'ytick',linspace(0,360,5),'FontName','times new roman','FontSize',12)
                end
                
                hold on              
                
            end            
            
            
            if (LOG_OBS_VEL=='T') %����������֤���ͼ���أ�T����F������
                
                TSR_Obs_Filepath='./Observation_Data/Velocity_Sequence'; %վλ��֤����·����·���ֶ�����
                TSR_Obs_Filename=[Sitename '.dat']; %վλ��֤���ļ���
                
                if exist(fullfile(TSR_Obs_Filepath,TSR_Obs_Filename)) %�������ƥ���վλ��֤�������򻭳���֤�㣬û�в���
                    %                     eval(['fid_TSR_Obs = fopen(''' fullfile(TSR_Obs_Filepath,TSR_Obs_Filename) ''');']); %����֤�ļ�
                    %                     TSR_Obsdata = fscanf(fid_TSR_Obs,'%f %f %f %f %f %f %f %f %f %f %f %f %f',[13 inf]); %��ȡ��֤������
                    %                     fclose(fid_TSR_Obs); %�ر��ļ�
                    
                    TSR_Obsdata = load(fullfile(TSR_Obs_Filepath,TSR_Obs_Filename));
                    
                    for i=1:Num_TSR_LAYER_VEL %������Ҫ��ͼ�Ĳ�������ѭ��
                        
                        eval(['subplot(' int2str(Num_TSR_LAYER_VEL*2) '1' int2str((i-1)*2+1) '),']);
                        %                         plot((TSR_Obsdata(1,:)/24),TSR_Obsdata(1+(TSR_LAYER_VEL(i)-1)*2+1,:),'b') %ʵ���i����������
                        %                         plot((TSR_Obsdata(1,:)/24),TSR_Obsdata(1+(TSR_LAYER_VEL(i)-1)*2+3,:),'r') %ʵ���i+1����������
                        plot((TSR_Obsdata(:,1)/24),TSR_Obsdata(:,1+(TSR_LAYER_VEL(i)-1)*2+1),'b') %ʵ���i����������
                        plot((TSR_Obsdata(:,1)/24),TSR_Obsdata(:,1+(TSR_LAYER_VEL(i)-1)*2+3),'r') %ʵ���i+1����������
                        
                        if(i==1) %��һ��ͼ����title
                            eval(['title(''' Sitename ' Velocity Validation'');'])    
                        end
                        
                        hold on
                        
                        eval(['subplot(' int2str(Num_TSR_LAYER_VEL*2) '1' int2str((i-1)*2+2) '),']);
                        %                         plot((TSR_Obsdata(1,:)/24),TSR_Obsdata(1+(TSR_LAYER_VEL(i)-1)*2+2,:),'b') %ʵ���i����������
                        %                         plot((TSR_Obsdata(1,:)/24),TSR_Obsdata(1+(TSR_LAYER_VEL(i)-1)*2+4,:),'r') %ʵ���i+1����������
                        plot((TSR_Obsdata(:,1)/24),TSR_Obsdata(:,1+(TSR_LAYER_VEL(i)-1)*2+2),'b') %ʵ���i����������
                        plot((TSR_Obsdata(:,1)/24),TSR_Obsdata(:,1+(TSR_LAYER_VEL(i)-1)*2+4),'r') %ʵ���i+1����������
                        
                        
                        hold on              
                    end            
                    
                end
                
            end
            
            eval(['print(gcf,''-dpng'',''' fullfile(TSR_Pic_Path,['Vel_' Sitename]) ''');'])
            
        end
        %---------------------------End Veloci Drawing-------------------------------           
        
        
        
        
        
        
        %---------------------------Salinty Drawing----------------------------------    
        if (LOG_TSR_S=='T') %�ж��Ƿ��ζȹ�����
            close(figure(1))
            disp('  Salinity Timeseries Drawing ...')
            disp(['    LAYER   ' num2str(TSR_LAYER_S)])
            Num_TSR_LAYER_S=length(TSR_LAYER_S); %������Ҫ������
            posf = [350 36 600 900];
            set(gcf,'position',[posf(1) posf(2) posf(3) posf(4)]); %���û�����λ�úʹ�С
            set(gcf,'paperpositionmode','auto');  %���ջ�����С����ͼ��
            set(gcf,'inverthardcopy','on'); %����ͼ��ʱȥ����ͼ���ڵĻ�ɫ
            
            for i=1:Num_TSR_LAYER_S %������Ҫ��ͼ�Ĳ�������ѭ��
                eval(['subplot(' int2str(Num_TSR_LAYER_S) '1' int2str(i) '),']); 
                plot((TSR_Compdata(1,:)/24)-TSR_LAG,TSR_Compdata(2+(TSR_LAYER_S(i)-1)*3+3,:),'k') %���Ƽ�����ζȵ�i�������ͼ
                set(gca,'FontName','times new roman','FontSize',12)
                if(i==1) %��һ��ͼ����title
                    eval(['title(''' Sitename ' Salinity Time Series'');'])    
                end
                
                
                if(i==Num_TSR_LAYER_S )  %���һ��ͼ����x����Ϣ                 
                    set(gca,'xlim',[TSR_BEG TSR_END],'xtick',linspace(TSR_BEG,TSR_END,TSR_END-TSR_BEG+1),'FontName','times new roman','FontSize',12)
                    xlabel('Time (Day)');        
                else
                    
                    set(gca,'xlim',[TSR_BEG TSR_END],'xtick',linspace(TSR_BEG,TSR_END,TSR_END-TSR_BEG+1),'xticklabel',[],'FontName','times new roman','FontSize',12)
                end
                ylabel('Salinty');
                
                hold on
                
            end
            
            if (LOG_OBS_S=='T') %վλ��֤���ͼ���أ�T����F������
                
                TSR_Obs_Filepath='./Observation_Data/Salinity_Sequence'; %վλ��֤����·����·���ֶ�����
                TSR_Obs_Filename=[Sitename '.dat']; %վλ��֤���ļ���
                
                if exist(fullfile(TSR_Obs_Filepath,TSR_Obs_Filename)) %�������ƥ���վλ��֤�������򻭳���֤�㣬û�в���
                    %                     eval(['fid_TSR_Obs = fopen(''' fullfile(TSR_Obs_Filepath,TSR_Obs_Filename) ''');']); %����֤�ļ�
                    %                     TSR_Obsdata = fscanf(fid_TSR_Obs,'%f %f %f %f %f %f %f',[7 inf]); %��ȡ��֤������
                    %                     fclose(fid_TSR_Obs); %�ر��ļ�
                    
                    TSR_Obsdata = load(fullfile(TSR_Obs_Filepath,TSR_Obs_Filename));
                    
                    for i=1:Num_TSR_LAYER_S %������Ҫ��ͼ�Ĳ�������ѭ��
                        
                        eval(['subplot(' int2str(Num_TSR_LAYER_S) '1' int2str((i)) '),']);
                        %                        plot((TSR_Obsdata(1,:)/24),TSR_Obsdata(1+TSR_LAYER_S(i),:),'b') %ʵ���i���ζ�����
                        %                        plot((TSR_Obsdata(1,:)/24),TSR_Obsdata(1+TSR_LAYER_S(i)+1,:),'r') %ʵ���i+1���ζ�����
                        plot((TSR_Obsdata(:,1)/24),TSR_Obsdata(:,1+TSR_LAYER_S(i)),'b') %ʵ���i���ζ�����
                        plot((TSR_Obsdata(:,1)/24),TSR_Obsdata(:,1+TSR_LAYER_S(i)+1),'r') %ʵ���i+1���ζ�����
                        set(gca,'FontName','times new roman','FontSize',12)
                        
                        if(i==1) %��һ��ͼ����title
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

if (LOG_TSR_SEC=='T') %����ʱ�����л�ͼ���أ�T����F������
    
    disp('================== Section Timeseries Drawing ====================')
    TSR_SEC_Pic_Path=fullfile(OUT_DIRE,'secflux'); %��ͼ������·��
    [s,mess,messid]=mkdir(TSR_SEC_Pic_Path); %������ͼ�������ļ���
    TSR_SEC_Comp_Filepath=fullfile(IN_DIRE,'secflux'); %ģʽ�������ļ�·��
    TSR_SEC_Comp_Fileinfo=dir(fullfile(TSR_SEC_Comp_Filepath,'*.out')); %ģʽ�������ļ���Ϣ�������ļ������޸�ʱ�䣬��С���Ƿ�Ϊ�ļ�·����
    TSR_SEC_Num=length(TSR_SEC_Comp_Fileinfo); %�����м�������ͼ���ļ�
    
    for i=1:TSR_SEC_Num %�����ļ���������ѭ����ͼ
        
        close(figure(1))
        TSR_SEC_Comp_Filename=TSR_SEC_Comp_Fileinfo(i).name; %�������ļ����ļ���
        Secname=TSR_SEC_Comp_Filename(9:length(TSR_SEC_Comp_Filename)-4); %���ļ�����ȡ��վ������
        disp(['  Drawing Section ' Secname ' ...'])
        eval(['fid_TSR_SEC_Comp = fopen(''' fullfile(TSR_SEC_Comp_Filepath,TSR_SEC_Comp_Filename) ''');']); %�򿪼������ļ�
        
        for i=1:4 %���ļ���ע����Ϣ
            NOTE=fgetl(fid_TSR_SEC_Comp);
        end
        
        TSR_SEC_Compdata = fscanf(fid_TSR_SEC_Comp,'%f %f %f %f %f',[5 inf]); %��������
        fclose(fid_TSR_SEC_Comp); %�ر�ģʽ�������ļ�
        
        posf = [350 36 600 900];
        set(gcf,'position',[posf(1) posf(2) posf(3) posf(4)]); %���û�����λ�úʹ�С
        set(gcf,'paperpositionmode','auto');  %���ջ�����С����ͼ��
        set(gcf,'inverthardcopy','on'); %����ͼ��ʱȥ����ͼ���ڵĻ�ɫ
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

