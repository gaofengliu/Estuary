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

    Init_day=datenum([IYEAR IMONTH IDAY0]); %����ģʽ���������    
    %----------��ȡ�������������ˮ������-------------------------
    fid_XYH_BIN=fopen('ch_hzbc_griddepth'); %�������������ꡢˮ���ļ������ļ���ģʽ���
    brecord=fread(fid_XYH_BIN,1,'integer*4'); %����¼��С��Ϣ
    h=fread(fid_XYH_BIN,IM*JM,'real*4'); %��ˮ������
    xr=fread(fid_XYH_BIN,IM*JM,'real*4'); %���������ĵ��X����
    yr=fread(fid_XYH_BIN,IM*JM,'real*4'); %���������ĵ��Y����
    erecord=fread(fid_XYH_BIN,1,'integer*4'); %����¼��С��Ϣ
    fclose(fid_XYH_BIN);
    
    if (brecord~=erecord) %���һ����¼��ǰ������¼�Ĵ�С��һ�£����˳�����
        disp('Error in reading file ch_hzbc_griddepth, please check!')
%       break
    end
    
    
    h=reshape(h,IM,JM); %ת������Ӧ�ľ�����ʽ
    xr=reshape(xr,IM,JM); %ת������Ӧ�ľ�����ʽ
    yr=reshape(yr,IM,JM); %ת������Ӧ�ľ�����ʽ
    
    
    fid_XYH_XY=fopen('ch_hzbc_xy_h.txt','w'); %�����ļ�
    for j=1:JM
        for i=1:IM
            fprintf(fid_XYH_XY,'%20.6f %19.6f %19.6f\r\n',xr(i,j),yr(i,j),h(i,j)); %����xr(i,j),yr(i,j),h(i,j)�ĸ�ʽд��������Ϣ
        end
    end
    fclose(fid_XYH_XY);
    
    if (FPT_COORDINATE=='XY') %���ƽ��ͼ��Ϊ54����
        xr=xr/1000; %��54���굥λ��mת��Ϊkm
        yr=yr/1000; %��54���굥λ��mת��Ϊkm
    end
    
    if (FPT_COORDINATE=='BL') %���ƽ��ͼ����Ҫ�Ծ�γ��Ϊ����
        disp('Please Trans  ''ch_hzbc_xy_h.txt''  To  ''ch_hzbc_xy_h.dat''  With  ''Coor.exe''') %������ת�����߽�������ת��
        pause
        
        
        fid_XYH_BL=fopen('ch_hzbc_xy_h.dat','r');
        XYH_BL=fscanf(fid_XYH_BL,'%f %f %f',[3 inf]); %����ת���õľ�γ�������ˮ������
        xr=reshape(XYH_BL(1,:),IM,JM); %������д�ɾ�����ʽ
        yr=reshape(XYH_BL(2,:),IM,JM); %��γ��д�ɾ�����ʽ       
        fclose(fid_XYH_BL);
    end
    %----------------------END------------------------
    
    
    %----------��ȡ��ͼ����İ�������-------------------------
    nisland=13;
    eval(['load LandFiles_' FPT_COORDINATE  '\land.dat']);
    for k1 = 1:nisland;
        eval(['load LandFiles_' FPT_COORDINATE '\island' int2str(k1) '.dat']);
    end
    
    eval(['load LandFiles_' FPT_COORDINATE  '\nanhuibiantan_weiken.dat']);
    
    eval(['load LandFiles_' FPT_COORDINATE  '\deep_water_way.dat']);
    
    if (FPT_COORDINATE=='XY') %���ʹ��54���꣬���������굥λ��mת��Ϊkm
        land=land/1000;
        
        for k1 = 1:nisland;
            eval(['island' int2str(k1) '=island' int2str(k1) '/1000;']);
        end
        
        nanhuibiantan_weiken=nanhuibiantan_weiken/1000;
        deep_water_way=deep_water_way/1000;
    end
    
    %-----------�·ݵ�Ӣ�ļ�д---------------------
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
    
    
    
    if (FPT_COLOR_TYPE=='C') %���ջ�ͼ����ɫ�������õ�ͼ�ı���ɫ
        bgcolor=[0 1 0];
    elseif (FPT_COLOR_TYPE=='G')
        bgcolor=[0.3 0.3 0.3];
    end
    
    %------------- Output Site and Section Position Drawing ---------------
    if (LOG_OPT=='T') %�����Ҫ�����վ��Ͷ���ʾ��ͼ
        disp('Output Site and Section Position Drawing ...')
        close(figure(1))
        OPT_Pic_Path=fullfile(OUT_DIRE,'output'); %���վ��Ͷ���ʾ��ͼ���·��
        [s,mess,messid]=mkdir(OPT_Pic_Path); %���վ��Ͷ���ʾ��ͼ��ŵ��ļ���
        OPT_Comp_Site_Filepath=fullfile(IN_DIRE,'timeseries'); %ģʽ��վ��������ļ�·��
        OPT_Comp_Section_Filepath=fullfile(IN_DIRE,'secflux'); %ģʽ�Ķ���������ļ�·��
        
        OPT_Comp_Site_Fileinfo=dir(fullfile(OPT_Comp_Site_Filepath,'*.out')); %ģʽ����վ���ļ���Ϣ�������ļ������޸�ʱ�䣬��С���Ƿ�Ϊ�ļ�·����
        OPT_Comp_Section_Fileinfo=dir(fullfile(OPT_Comp_Section_Filepath,'*.out')); %ģʽ��������ļ���Ϣ�������ļ������޸�ʱ�䣬��С���Ƿ�Ϊ�ļ�·����
        OPT_Site_Num=length(OPT_Comp_Site_Fileinfo); %�����м��������        
        OPT_Section_Num=length(OPT_Comp_Section_Fileinfo); %�����м����������    
        
      for i=1:OPT_Section_Num
            OPT_Comp_Section_Filename=OPT_Comp_Section_Fileinfo(i).name; %�������ļ����ļ���
            eval(['fid_OPT_Comp_Section = fopen(''' fullfile(OPT_Comp_Section_Filepath,OPT_Comp_Section_Filename) ''');']); %��վ�����ļ�
            
            
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
            
            fclose(fid_OPT_Comp_Section); %�ر�ģʽ�������ļ�
            
            plot([xr(i1,j1) xr(i2,j2)],[yr(i1,j1),yr(i2,j2)],'linewidth',1.5,'color','k');
            
            hold on
            text(xr(i1,j1),yr(i1,j1),secname(5:end),'FontName','times new roman','FontSize',12)
            
        end
        
        
        
        
        set(gca,'box','on','FontName','times new roman','FontSize',12);
        
        set(gca,'color',[1 1 1]);
        
        set(gca,'xlim',[FPT_XMIN FPT_XMAX],'xtick',linspace(FPT_XMIN,FPT_XMAX,4));
        
        set(gca,'ylim',[FPT_YMIN FPT_YMAX],'ytick',linspace(FPT_YMIN,FPT_YMAX,4));
        
        set(gcf,'inverthardcopy','off'); %����ͼƬʱ�������õ���ɫ�����Զ�����
        
        set(gcf,'color',[1 1 1]);
        
        hold on 
     
        %---------���ư��ߡ���------------------------
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
        
        
        if (FPT_COORDINATE=='BL') %����ͼ���������ѡȡ������Ӧ��label
            xlabel('Longitude (\circE)');
            ylabel('Latitude (\circN)');
        elseif (FPT_COORDINATE=='XY')
            xlabel('Distance (km)');
            ylabel('Distance (km)');
        end     
        
        eval(['print(gcf,''-depsc'',''' fullfile(OPT_Pic_Path,'\Output') ''');']); %����ͼ��
        
    end
    
    
    
    
    
    %------------- Elevation Field Distribution Drawing -------------------
    if (LOG_EPT=='T') %�����Ҫ��ˮλ��ƽ��ͼ
        disp('Elevation Field Distribution Drawing ...')
        EPT_Pic_Path=fullfile(OUT_DIRE,'field_distri\elevation'); %ˮλƽ��ͼ���·��
        [s,mess,messid]=mkdir(EPT_Pic_Path); %����ˮλƽ��ͼ��ŵ��ļ���
        EPT_Comp_Filepath=fullfile(IN_DIRE,'field_distri\elevation'); %ģʽ��ˮλ���������ļ�·��
        
        
        EPT_TIME=str2num(EPT_TIME); %�������ˮλ��ͼʱ�̵��ַ�����ת��Ϊ��������
        Num_EPT_TIME=length(EPT_TIME); %�ܹ���Ҫ������ʱ�̵�ͼ��
        EPT_contourscale=linspace(-5,5,41); %��Ҫ���ĵ�ֵ��
        
        for i=1:Num_EPT_TIME %������Ҫ��ͼ����������ѭ��
            close(figure(1));
            disp(['  The ' num2str(EPT_TIME(i)) ' Hours Drawing'])
            %------------------ Elevation Field Distribution Data Reading ------------------------
            EPT_Comp_Filename_SN=EPT_TIME(i)/(N_FPT/3600); %����ģʽ�����ʱ������������Ҫ��ͼʱ�̶�Ӧ�����
            EPT_Comp_Filename=sprintf('el_field_%06.6d',EPT_Comp_Filename_SN); %ģʽ���������ˮλ�ļ�
            
            if ~exist(fullfile(EPT_Comp_Filepath,EPT_Comp_Filename)) %�ж�ˮλ���ļ��Ƿ���ڣ��������������ʾ��Ϣ
                disp(['    ' EPT_Comp_Filename ' dose not exist, drawing failed, please check the setting!']) %�����ļ������ڣ������ʾ��Ϣ
                
            else %��������ļ���������л�ͼ
                eval(['fid_EPT_Comp = fopen(''' fullfile(EPT_Comp_Filepath,EPT_Comp_Filename) ''',''r'',''b'');']); %���ļ������ڸö������ļ����á�BIG_ENDIAN����װ�����Դ���ѡ�������b��
                brecord=fread(fid_EPT_Comp,1,'integer*4'); %����¼��С��Ϣ
                thour=fread(fid_EPT_Comp,1,'real*4');
                EPT_Comp_Data=fread(fid_EPT_Comp,IM*JM,'real*4');
                erecord=fread(fid_EPT_Comp,1,'integer*4'); %����¼��С��Ϣ
                
                if (brecord~=erecord) %���һ����¼��ǰ������¼�Ĵ�С��һ�£����˳�����
                    disp(['Error in reading file ' fullfile(EPT_Comp_Filepath,EPT_Comp_Filename) ', please check!'])
                    break
                end
                EPT_Comp_Data=reshape(EPT_Comp_Data,IM,JM); %��ˮλ����ת��Ϊ������ʽ
                fclose(fid_EPT_Comp); %�رռ��������ˮλ�ļ�
                %---------------- End Elevation Field Distribution Data Reading ----------------------
                
                for i1=1:IM
                    for j1=1:JM
                        if (h(i1,j1)+EPT_Comp_Data(i1,j1)<DMIN||h(i1,j1)<-10) %����̲��½�ش���ˮλ��ֵΪnan
                            EPT_Comp_Data(i1,j1)=NaN;
                        end
                    end
                end  
                
                if (FPT_COLOR_TYPE=='C') %���ݻ�ͼ����ɫ���ͣ�ѡ����ɫ���ļ�
                    map=load ('Colormap\cm_el_C.dat');
                elseif (FPT_COLOR_TYPE=='G')
                    map=load ('Colormap\cm_el_G.dat');
                end
                map=[bgcolor;map]; %��ӱ�����ɫ����ˮλΪnan��λ�ð�����ɫ��
                colormap(map);  %������ɫ��
                
                Pic_time=datevec(Init_day+thour/24); %����ģʽ���õ���ʼʱ��͸����ݼ�¼��thour����ʵ�ʵ�ʱ�䣬Pic_time���� �ꡢ�¡��ա�ʱ���֡��� ��Ϣ
                
                set(gca,'box','on','FontName','times new roman','FontSize',12);
                set(gca,'color',bgcolor);
                set(gca,'xlim',[FPT_XMIN FPT_XMAX],'xtick',linspace(FPT_XMIN,FPT_XMAX,4));
                set(gca,'ylim',[FPT_YMIN FPT_YMAX],'ytick',linspace(FPT_YMIN,FPT_YMAX,4));
                set(gcf,'inverthardcopy','off'); %����ͼƬʱ�������õ���ɫ�����Զ�����
                set(gcf,'color',[1 1 1]);
                hold on 
                
                warning off; %��contourfʱ���кܶྯ����Ϣ������Ӱ�컭ͼ�����Թص�warning
                contourf(xr,yr,EPT_Comp_Data,EPT_contourscale);
                shading flat; 
                
                EPT_Comp_Data_Max=max(max(EPT_Comp_Data)); %���ˮλ���ݵ����ֵ��Сֵ
                EPT_Comp_Data_Min=min(min(EPT_Comp_Data));
                cmax=ceil(max(abs(EPT_Comp_Data_Max),abs(EPT_Comp_Data_Min))); %������������������ȡ��
                cmin=-cmax; %��Ϊ��ɫ�������ɫ�ǶԳƵģ�����Ҫʹ��������Գ�
                caxis([cmin-(cmax-cmin)/47 cmax]); %Ϊ�˰ѱ���ɫ�ŵ���ɫ�������ֲ���ʾ��colorbar���47����������ɫ������47����ɫ
                
                if (FPT_COLOR_TYPE=='C') %�������ɫͼ������Ҫ��colorbar
                    hc=colorbar;
                    set(hc,'ylim',[cmin cmax],'FontName','times new roman','FontSize',12); %����colorbar����ʾ��Χ
                end
                
                [cc,hh]=contour(xr,yr,EPT_Comp_Data,EPT_contourscale);
                clabel(cc,hh,'FontSize',10,'FontName','Times new roman');
                set(hh,'EdgeColor',[0.3 0.3 0.3],'linewidth',1) %���õ�ֵ�ߵ���ɫ�ʹ�ϸ
                
                % 				title([num2str(Pic_time(1)) '��' num2str(Pic_time(2)) '��' num2str(Pic_time(3)) '��' num2str(Pic_time(4)) 'ʱ ˮλ�ֲ�'])
                title([num2str(Pic_time(4)) ':00 ' mon(Pic_time(2),:) ' ' num2str(Pic_time(3)) ' ' num2str(Pic_time(1)) ' (GMT +8)'])    
                
                hold on  
                
                %---------���ư��ߡ���------------------------
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
                
                if (FPT_COORDINATE=='BL') %����ͼ���������ѡȡ������Ӧ��label
                    xlabel('Longitude (\circE)');
                    ylabel('Latitude (\circN)');
                elseif (FPT_COORDINATE=='XY')
                    xlabel('Distance (km)');
                    ylabel('Distance (km)');
                end
                
                eval(['print(gcf,''-dpng'',''' fullfile(EPT_Pic_Path,EPT_Comp_Filename) ''');']); %����ͼ��
                
            end
            
        end      
    end    
    %------------- End Elevation Field Distribution Drawing --------------
    
    %------------- Velocity Field Distribution Drawing -------------------
    if (LOG_VPT_UV=='T'||LOG_VPT_SD=='T') %�����Ҫ��������ƽ��ͼ
        disp('Velocity Field Distribution Drawing ...')
        VPT_Pic_Path=fullfile(OUT_DIRE,'field_distri\current'); %�ζ�ƽ��ͼ���·��
        [s,mess,messid]=mkdir(VPT_Pic_Path); %�����ζ�ƽ��ͼ�Ĵ���ļ���
        EPT_Comp_Filepath=fullfile(IN_DIRE,'field_distri\elevation'); %ģʽ��ˮλ���������ļ�·��
        VPT_Comp_Filepath=fullfile(IN_DIRE,'field_distri\current'); %ģʽ���ζȳ��������ļ�·��
        
        
        VPT_TIME=str2num(VPT_TIME); %��������ζȻ�ͼʱ�̵��ַ�����ת��Ϊ��������
        Num_VPT_TIME=length(VPT_TIME); %�ܹ���Ҫ������ʱ�̵�ͼ��
        VPT_LAYER=str2num(VPT_LAYER); %��������ζȻ�ͼ�������ַ�����ת��Ϊ��������
        Num_VPT_LAYER=length(VPT_LAYER); %�ܹ���Ҫ��������ε�ͼ��
        H_contourscale=[-1 0 5 10 15 20 25 30 50 70 100]; %ʸ����ͼ��ˮ���ֵ������
        SD_contourscale=[-1 0 5 10 15 20 25 30 50 70 100]; %������ֵ��ֵ������
        
        
        
        for i=1:Num_VPT_TIME %������Ҫ��ͼ����������ѭ��
            
            disp(['  The ' num2str(VPT_TIME(i)) ' Hours Drawing'])
            
            %------ Elevation Field Distribution Data Reading ,For Check Tideflat Grid ----------
            EPT_Comp_Filename_SN=VPT_TIME(i)/(N_FPT/3600); %����ģʽ�����ʱ������������Ҫ��ͼʱ�̶�Ӧ�����
            EPT_Comp_Filename=sprintf('el_field_%06.6d',EPT_Comp_Filename_SN); %ģʽ���������ˮλ�ļ�
            VPT_Comp_Filename_SN=VPT_TIME(i)/(N_FPT/3600); %����ģʽ�����ʱ������������Ҫ��ͼʱ�̶�Ӧ�����
            VPT_Comp_Filename=sprintf('v_field_%06.6d',VPT_Comp_Filename_SN); %ģʽ���������ˮλ�ļ�           
            
            if ~(exist(fullfile(EPT_Comp_Filepath,EPT_Comp_Filename))&&exist(fullfile(VPT_Comp_Filepath,VPT_Comp_Filename))) %�ж�ˮλ���ļ������ٳ��ļ��Ƿ񶼴��ڣ��������������ʾ��Ϣ
                disp(['    ' EPT_Comp_Filename ' or ' VPT_Comp_Filename ' dose not exist, drawing failed, please check the setting!']) %�����ļ������ڣ������ʾ��Ϣ
            else %�����ļ����ڣ����ͼ
                eval(['fid_EPT_Comp = fopen(''' fullfile(EPT_Comp_Filepath,EPT_Comp_Filename) ''',''r'',''b'');']); %���ļ������ڸö������ļ����á�BIG_ENDIAN����װ�����Դ���ѡ�������b��
                brecord=fread(fid_EPT_Comp,1,'integer*4'); %����¼��С��Ϣ
                thour=fread(fid_EPT_Comp,1,'real*4');
                EPT_Comp_Data=fread(fid_EPT_Comp,IM*JM,'real*4');
                erecord=fread(fid_EPT_Comp,1,'integer*4'); %����¼��С��Ϣ
                
                if (brecord~=erecord) %���һ����¼��ǰ������¼�Ĵ�С��һ�£����˳�����
                    disp(['Error in reading file ' fullfile(EPT_Comp_Filepath,EPT_Comp_Filename) ', please check!'])
                    break
                end
                EPT_Comp_Data=reshape(EPT_Comp_Data,IM,JM); %��ˮλ����ת��Ϊ������ʽ
                fclose(fid_EPT_Comp); %�رռ��������ˮλ�ļ�
                %---------------- End Elevation Field Distribution Data Reading ----------------------
                
                
                
                %--------------------- Velocity Field Distribution Data Reading -----------------------
                eval(['fid_VPT_Comp = fopen(''' fullfile(VPT_Comp_Filepath,VPT_Comp_Filename) ''',''r'',''b'');']); %���ļ������ڸö������ļ����á�BIG_ENDIAN����װ�����Դ���ѡ�������b��
                brecord=fread(fid_VPT_Comp,1,'integer*4'); %����¼��С��Ϣ
                thour=fread(fid_VPT_Comp,1,'real*4');
                VPT_Comp_U=fread(fid_VPT_Comp,IM*JM*KB,'real*4');
                VPT_Comp_V=fread(fid_VPT_Comp,IM*JM*KB,'real*4');
                erecord=fread(fid_VPT_Comp,1,'integer*4'); %����¼��С��Ϣ
                
                if (brecord~=erecord) %���һ����¼��ǰ������¼�Ĵ�С��һ�£����˳�����
                    disp(['Error in reading file ' fullfile(VPT_Comp_Filepath,VPT_Comp_Filename) ', please check!'])
                    break
                end
                VPT_Comp_U=reshape(VPT_Comp_U,IM,JM,KB); %������Uת��Ϊ������ʽ
                VPT_Comp_V=reshape(VPT_Comp_V,IM,JM,KB); %������Vת��Ϊ������ʽ
                fclose(fid_VPT_Comp); %�رռ�������������ļ�
                %------------------- End Velocity Field Distribution Data Reading --------------------
                
                
                
                %��ʼ�����ڻ�������ʸ��ͼ�ı���������ֵΪnan
                X_VPT=NaN*ones(IM,JM);
                Y_VPT=NaN*ones(IM,JM);
                U_VPT=NaN*ones(IM,JM,KB);
                V_VPT=NaN*ones(IM,JM,KB);
                
                %�������õĻ�ͼ���VPT_INTERVAL������Ч���λ����������ٸ�ֵ
                for i2=1:VPT_INTERVAL:IM 
                    for j2=1:VPT_INTERVAL:JM
                        X_VPT(i2,j2)=xr(i2,j2);
                        Y_VPT(i2,j2)=yr(i2,j2);                       
                        U_VPT(i2,j2,:)=VPT_Comp_U(i2,j2,:);
                        V_VPT(i2,j2,:)=VPT_Comp_V(i2,j2,:);
                    end
                end    				
                
                
                
                H_VPT=h; %��ˮ�����ݸ�ֵ��һ����ʱ����H_VPT�����Խ�H_VPT���в���
                for i1=1:IM
                    for j1=1:JM
                        if (h(i1,j1)+EPT_Comp_Data(i1,j1)<DMIN||h(i1,j1)<-10) %����̲��½�ش���ˮ������ٸ�ֵΪnan
                            H_VPT(i1,j1)=NaN;
                            X_VPT(i1,j1)=NaN;
                            Y_VPT(i1,j1)=NaN;                    
                            U_VPT(i1,j1,:)=NaN; %���ڻ�ʸ��ͼ
                            V_VPT(i1,j1,:)=NaN;  
                            VPT_Comp_U(i1,j1,:)=NaN; %���ڻ���ֵͼ
                            VPT_Comp_V(i1,j1,:)=NaN;  
                            
                        end
                    end
                end  
                
                %Ϊ�˻�����ʸ��С��ǣ�������ĩβ���һ��ˮƽ���������
                X_VPT(i1,j1+1)=FPT_XMIN+(FPT_XMAX-FPT_XMIN)*0.05; %����Ϊ��������5%��λ��
                Y_VPT(i1,j1+1)=FPT_YMIN+(FPT_YMAX-FPT_YMIN)*0.1; %����Ϊ��С����10%��λ��
                U_VPT(i1,j1+1,:)=VPT_SCALE;
                V_VPT(i1,j1+1,:)=0;
                
                %Ϊ�˻�����ʸ��С��ǣ�������ĩβ���һ����ֱ���������
                X_VPT(i1+1,j1+1)=FPT_XMIN+(FPT_XMAX-FPT_XMIN)*0.05; %����Ϊ��������5%��λ��
                Y_VPT(i1+1,j1+1)=FPT_YMIN+(FPT_YMAX-FPT_YMIN)*0.12; %����Ϊ��С����12%��λ��
                U_VPT(i1+1,j1+1,:)=0;
                V_VPT(i1+1,j1+1,:)=VPT_SCALE;		
                
                %��quiver�������趨��ͷ��ʵ�ʴ�С�����������Զ����ڴ�С��������Ҫ�ѽ�ԭʼ���ݷŴ󵽺���������Ӧ�Ĵ�С�߶�
                U_VPT=U_VPT*(FPT_XMAX-FPT_XMIN)*0.05/VPT_SCALE; %(FPT_XMAX-FPT_XMIN)*0.05/VPT_SCALE��ʾ��x���5%�ĳ��ȶ�ӦVPT_SCALE�Ĵ�С
                V_VPT=V_VPT*(FPT_YMAX-FPT_YMIN)*0.05/VPT_SCALE; %(FPT_YMAX-FPT_YMIN)*0.05/VPT_SCALE��ʾ��y���5%�ĳ��ȶ�ӦVPT_SCALE�Ĵ�С
                
                Pic_time=datevec(Init_day+thour/24); %����ģʽ���õ���ʼʱ��͸����ݼ�¼��thour����ʵ�ʵ�ʱ�䣬Pic_time���� �ꡢ�¡��ա�ʱ���֡��� ��Ϣ
                
                
                for k=1:Num_VPT_LAYER
                    
                    if (LOG_VPT_UV=='T') %��������ƽ��ʸ��ͼ
                        
                        close(figure(1));
                        disp(['    Layer ' num2str(VPT_LAYER(k))])
                        if (FPT_COLOR_TYPE=='C') %���ݻ�ͼ����ɫ���ͣ�ѡ����ɫ���ļ�
                            map=load ('Colormap\cm_v_C.dat');
                        elseif (FPT_COLOR_TYPE=='G')
                            map=load ('Colormap\cm_v_G.dat');
                        end
                        map=[bgcolor;map]; %��ӱ�����ɫ�����ζ�Ϊnan��λ�ð�����ɫ��
                        colormap(map);  %������ɫ��
                        
                        
                        set(gca,'box','on','FontName','times new roman','FontSize',12);
                        set(gca,'color',bgcolor);
                        set(gca,'xlim',[FPT_XMIN FPT_XMAX],'xtick',linspace(FPT_XMIN,FPT_XMAX,4));
                        set(gca,'ylim',[FPT_YMIN FPT_YMAX],'ytick',linspace(FPT_YMIN,FPT_YMAX,4));
                        set(gcf,'inverthardcopy','off'); %����ͼƬʱ�������õ���ɫ�����Զ�����
                        set(gcf,'color',[1 1 1]);
                        hold on 
                        
                        
                        warning off; %��contourfʱ���кܶྯ����Ϣ������Ӱ�컭ͼ�����Թص�warning
                        contourf(xr,yr,H_VPT,H_contourscale);
                        shading flat; 
                        caxis([-1 70]); 
                        if (FPT_COLOR_TYPE=='C') %�������ɫͼ������Ҫ��colorbar
                            hc=colorbar;
                            set(hc,'ylim',[0 70],'FontName','times new roman','FontSize',12); %����colorbar����ʾ��Χ
                            [cc,hh]=contour(xr,yr,H_VPT,H_contourscale); %��ͼ��ˮ��ĵ�ֵ�ߣ��Ҷ�ͼ����
                            clabel(cc,hh,'FontSize',10,'FontName','Times new roman');
                            set(hh,'EdgeColor',[0.3 0.3 0.3],'linewidth',1) %���õ�ֵ�ߵ���ɫ�ʹ�ϸ
                        end
                        
                        eval(['quiver(X_VPT,Y_VPT,U_VPT(:,:,' num2str(VPT_LAYER(k)) '),V_VPT(:,:,' num2str(VPT_LAYER(k)) '),0,''k'');']);
                        
                        % 					if (VPT_LAYER(k)==1)
                        % 						VPT_LAYER_NAME='��';
                        % 					elseif (VPT_LAYER(k)==5)
                        % 						VPT_LAYER_NAME='��';
                        % 					else
                        % 						VPT_LAYER_NAME=['��' num2str(VPT_LAYER(k))];
                        % 					end
                        % 					title([num2str(Pic_time(1)) '��' num2str(Pic_time(2)) '��' num2str(Pic_time(3)) '��' num2str(Pic_time(4)) 'ʱ ' VPT_LAYER_NAME '�����ٷֲ�'])
                        
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
                        
                        %---------���ư��ߡ���------------------------
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
                        
                        if (FPT_COORDINATE=='BL') %����ͼ���������ѡȡ������Ӧ��label
                            xlabel('Longitude (\circE)');
                            ylabel('Latitude (\circN)');
                        elseif (FPT_COORDINATE=='XY')
                            xlabel('Distance (km)');
                            ylabel('Distance (km)');
                        end
                        
                        eval(['print(gcf,''-dpng'',''' fullfile(VPT_Pic_Path,VPT_Comp_Filename) '_UV_' num2str(VPT_LAYER(k)) ''');']); %����ͼ��
                        
                    end %��������ƽ��ʸ��ͼ end
                    
                    
                    
                    
                    if (LOG_VPT_SD=='T') %��������ƽ����ֵͼ����ֵ��ͼ��
                        
                        close(figure(1));
                        disp(['    Layer ' num2str(VPT_LAYER(k))])
                        if (FPT_COLOR_TYPE=='C') %���ݻ�ͼ����ɫ���ͣ�ѡ����ɫ���ļ�
                            map=load ('Colormap\cm_speed_C.dat');
                        elseif (FPT_COLOR_TYPE=='G')
                            map=load ('Colormap\cm_speed_G.dat');
                        end
                        map=[bgcolor;map]; %��ӱ�����ɫ�����ζ�Ϊnan��λ�ð�����ɫ��
                        colormap(map);  %������ɫ��
                        
                        
                        set(gca,'box','on','FontName','times new roman','FontSize',12);
                        set(gca,'color',bgcolor);
                        set(gca,'xlim',[FPT_XMIN FPT_XMAX],'xtick',linspace(FPT_XMIN,FPT_XMAX,4));
                        set(gca,'ylim',[FPT_YMIN FPT_YMAX],'ytick',linspace(FPT_YMIN,FPT_YMAX,4));
                        set(gcf,'inverthardcopy','off'); %����ͼƬʱ�������õ���ɫ�����Զ�����
                        set(gcf,'color',[1 1 1]);
                        hold on 
                        
                        VPT_Comp_Data(:,:,VPT_LAYER(k))=sqrt(VPT_Comp_U(:,:,VPT_LAYER(k)).^2+VPT_Comp_V(:,:,VPT_LAYER(k)).^2); %��uvת��Ϊ����
                        
                        warning off; %��contourfʱ���кܶྯ����Ϣ������Ӱ�컭ͼ�����Թص�warning
                        contourf(xr,yr,VPT_Comp_Data(:,:,VPT_LAYER(k)));
                        
                        shading flat; 
                        
                        
                        % 						caxis([-6/64 3]); %���������ֵΪ3�����ǰ������������������ֵָ�����ֵ
                        VPT_Comp_Data_Max=max(max(VPT_Comp_Data(:,:,VPT_LAYER(k)))); %����������ݵ����ֵ��Сֵ
                        VPT_Comp_Data_Min=min(min(VPT_Comp_Data(:,:,VPT_LAYER(k))));
                        cmax=ceil(max(abs(VPT_Comp_Data_Max),abs(VPT_Comp_Data_Min))); %������������������ȡ��
                        caxis([-cmax/64 cmax]); %Ϊ�˰ѱ���ɫ�ŵ���ɫ�������ֲ���ʾ��colorbar���64����������ɫ������64����ɫ                 
                        
                        if (FPT_COLOR_TYPE=='C') %�������ɫͼ������Ҫ��colorbar
                            hc=colorbar;
                            set(hc,'ylim',[0 cmax],'FontName','times new roman','FontSize',12); %����colorbar����ʾ��Χ
                            [cc,hh]=contour(xr,yr,VPT_Comp_Data(:,:,VPT_LAYER(k))); 
                            clabel(cc,hh,'FontSize',10,'FontName','Times new roman');
                            set(hh,'EdgeColor',[0.3 0.3 0.3],'linewidth',1) %���õ�ֵ�ߵ���ɫ�ʹ�ϸ
                        elseif (FPT_COLOR_TYPE=='G')
                            [cc,hh]=contour(xr,yr,VPT_Comp_Data(:,:,VPT_LAYER(k))); 
                            clabel(cc,hh,'FontSize',10,'FontName','Times new roman');
                            set(hh,'EdgeColor',[0.3 0.3 0.3],'linewidth',1) %���õ�ֵ�ߵ���ɫ�ʹ�ϸ
                        end
                        
                        % 						eval(['quiver(X_VPT,Y_VPT,U_VPT(:,:,' num2str(VPT_LAYER(k)) '),V_VPT(:,:,' num2str(VPT_LAYER(k)) '),0,''k'');']);
                        
                        % 					if (VPT_LAYER(k)==1)
                        % 						VPT_LAYER_NAME='��';
                        % 					elseif (VPT_LAYER(k)==5)
                        % 						VPT_LAYER_NAME='��';
                        % 					else
                        % 						VPT_LAYER_NAME=['��' num2str(VPT_LAYER(k))];
                        % 					end
                        % 					title([num2str(Pic_time(1)) '��' num2str(Pic_time(2)) '��' num2str(Pic_time(3)) '��' num2str(Pic_time(4)) 'ʱ ' VPT_LAYER_NAME '�����ٷֲ�'])
                        
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
                        
                        %---------���ư��ߡ���------------------------
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
                        
                        if (FPT_COORDINATE=='BL') %����ͼ���������ѡȡ������Ӧ��label
                            xlabel('Longitude (\circE)');
                            ylabel('Latitude (\circN)');
                        elseif (FPT_COORDINATE=='XY')
                            xlabel('Distance (km)');
                            ylabel('Distance (km)');
                        end
                        
                        eval(['print(gcf,''-dpng'',''' fullfile(VPT_Pic_Path,VPT_Comp_Filename) '_SD_' num2str(VPT_LAYER(k)) ''');']); %����ͼ��
                        
                    end %��������ƽ����ֵͼ����ֵ��ͼ��end
                    
                end
                clear H_VPT X_VPT Y_VPT U_VPT V_VPT
            end
        end
    end
    %------------- End Velocity Field Distribution Drawing ----------------    
    
    %------------- Salinity Field Distribution Drawing -------------------
    if (LOG_SPT=='T') %�����Ҫ���ζȵ�ƽ��ͼ
        disp('Salinity Field Distribution Drawing ...')
        SPT_Pic_Path=fullfile(OUT_DIRE,'field_distri\salinity'); %�ζ�ƽ��ͼ���·��
        [s,mess,messid]=mkdir(SPT_Pic_Path); %�����ζ�ƽ��ͼ�Ĵ���ļ���
        EPT_Comp_Filepath=fullfile(IN_DIRE,'field_distri\elevation'); %ģʽ��ˮλ���������ļ�·��
        SPT_Comp_Filepath=fullfile(IN_DIRE,'field_distri\salinity'); %ģʽ���ζȳ��������ļ�·��
        
        
        SPT_TIME=str2num(SPT_TIME); %��������ζȻ�ͼʱ�̵��ַ�����ת��Ϊ��������
        Num_SPT_TIME=length(SPT_TIME); %�ܹ���Ҫ������ʱ�̵�ͼ��
        SPT_LAYER=str2num(SPT_LAYER); %��������ζȻ�ͼ�������ַ�����ת��Ϊ��������
        Num_SPT_LAYER=length(SPT_LAYER); %�ܹ���Ҫ��������ε�ͼ��
        SPT_contourscale=[-0.5 0 0.5 1 3 5 10 15 20 25 30 35]; %��Ҫ���ĵ�ֵ��
        
        
        
        
        for i=1:Num_SPT_TIME %������Ҫ��ͼ����������ѭ��
            
            disp(['  The ' num2str(SPT_TIME(i)) ' Hours Drawing'])
            
            %------ Elevation Field Distribution Data Reading ,For Check Tideflat Grid ----------
            EPT_Comp_Filename_SN=SPT_TIME(i)/(N_FPT/3600); %����ģʽ�����ʱ������������Ҫ��ͼʱ�̶�Ӧ�����
            EPT_Comp_Filename=sprintf('el_field_%06.6d',EPT_Comp_Filename_SN); %ģʽ���������ˮλ�ļ�
            SPT_Comp_Filename_SN=SPT_TIME(i)/(N_FPT/3600); %����ģʽ�����ʱ������������Ҫ��ͼʱ�̶�Ӧ�����
            SPT_Comp_Filename=sprintf('s_field_%06.6d',SPT_Comp_Filename_SN); %ģʽ���������ˮλ�ļ�           
            
            if ~(exist(fullfile(EPT_Comp_Filepath,EPT_Comp_Filename))&&exist(fullfile(SPT_Comp_Filepath,SPT_Comp_Filename))) %�ж�ˮλ���ļ����ζȳ��ļ��Ƿ񶼴��ڣ��������������ʾ��Ϣ
                disp(['    ' EPT_Comp_Filename ' or ' SPT_Comp_Filename ' dose not exist, drawing failed, please check the setting!']) %�����ļ������ڣ������ʾ��Ϣ
            else %�����ļ����ڣ����ͼ
                eval(['fid_EPT_Comp = fopen(''' fullfile(EPT_Comp_Filepath,EPT_Comp_Filename) ''',''r'',''b'');']); %���ļ������ڸö������ļ����á�BIG_ENDIAN����װ�����Դ���ѡ�������b��
                brecord=fread(fid_EPT_Comp,1,'integer*4'); %����¼��С��Ϣ
                thour=fread(fid_EPT_Comp,1,'real*4');
                EPT_Comp_Data=fread(fid_EPT_Comp,IM*JM,'real*4');
                erecord=fread(fid_EPT_Comp,1,'integer*4'); %����¼��С��Ϣ
                
                if (brecord~=erecord) %���һ����¼��ǰ������¼�Ĵ�С��һ�£����˳�����
                    disp(['Error in reading file ' fullfile(EPT_Comp_Filepath,EPT_Comp_Filename) ', please check!'])
                    break
                end
                EPT_Comp_Data=reshape(EPT_Comp_Data,IM,JM); %��ˮλ����ת��Ϊ������ʽ
                fclose(fid_EPT_Comp); %�رռ��������ˮλ�ļ�
                %---------------- End Elevation Field Distribution Data Reading ----------------------
                
                
                
                %--------------------- Salinity Field Distribution Data Reading -----------------------
                eval(['fid_SPT_Comp = fopen(''' fullfile(SPT_Comp_Filepath,SPT_Comp_Filename) ''',''r'',''b'');']); %���ļ������ڸö������ļ����á�BIG_ENDIAN����װ�����Դ���ѡ�������b��
                brecord=fread(fid_SPT_Comp,1,'integer*4'); %����¼��С��Ϣ
                thour=fread(fid_SPT_Comp,1,'real*4');
                SPT_Comp_Data=fread(fid_SPT_Comp,IM*JM*KB,'real*4');
                erecord=fread(fid_SPT_Comp,1,'integer*4'); %����¼��С��Ϣ
                
                if (brecord~=erecord) %���һ����¼��ǰ������¼�Ĵ�С��һ�£����˳�����
                    disp(['Error in reading file ' fullfile(SPT_Comp_Filepath,SPT_Comp_Filename) ', please check!'])
                    break
                end
                SPT_Comp_Data=reshape(SPT_Comp_Data,IM,JM,KB); %���ζ�ת��Ϊ������ʽ
                fclose(fid_SPT_Comp); %�رռ��������ˮλ�ļ�
                %------------------- End Salinity Field Distribution Data Reading -----------------------
                
                
                for i1=1:IM
                    for j1=1:JM
                        if (h(i1,j1)+EPT_Comp_Data(i1,j1)<DMIN||h(i1,j1)<-10) %����̲��½�ش����ζȸ�ֵΪnan
                            SPT_Comp_Data(i1,j1,:)=NaN;
                        end
                        
                        for k=1:KB-1
                            if SPT_Comp_Data(i1,j1,k)<0
                                SPT_Comp_Data(i1,j1,k)=0;
                            end
                        end
                    end
                end  
                
                
                
                Pic_time=datevec(Init_day+thour/24); %����ģʽ���õ���ʼʱ��͸����ݼ�¼��thour����ʵ�ʵ�ʱ�䣬Pic_time���� �ꡢ�¡��ա�ʱ���֡��� ��Ϣ
                
                
                for k=1:Num_SPT_LAYER
                    close(figure(1));
                    disp(['    Layer ' num2str(SPT_LAYER(k))])
                    if (FPT_COLOR_TYPE=='C') %���ݻ�ͼ����ɫ���ͣ�ѡ����ɫ���ļ�
                        map=load ('Colormap\cm_s_C.dat');
                    elseif (FPT_COLOR_TYPE=='G')
                        map=load ('Colormap\cm_s_G.dat');
                    end
                    map=[bgcolor;map]; %��ӱ�����ɫ�����ζ�Ϊnan��λ�ð�����ɫ��
                    colormap(map);  %������ɫ��
                    
                    
                    set(gca,'box','on','FontName','times new roman','FontSize',12);
                    set(gca,'color',bgcolor);
                    set(gca,'xlim',[FPT_XMIN FPT_XMAX],'xtick',linspace(FPT_XMIN,FPT_XMAX,4));
                    set(gca,'ylim',[FPT_YMIN FPT_YMAX],'ytick',linspace(FPT_YMIN,FPT_YMAX,4));
                    set(gcf,'inverthardcopy','off'); %����ͼƬʱ�������õ���ɫ�����Զ�����
                    set(gcf,'color',[1 1 1]);
                    hold on 
                    
                    
                    warning off; %��contourfʱ���кܶྯ����Ϣ������Ӱ�컭ͼ�����Թص�warning
                    eval(['contourf(xr,yr,SPT_Comp_Data(:,:,' num2str(SPT_LAYER(k)) '),SPT_contourscale);']);
                    shading flat; 
                    caxis([-0.5 35]); %Ϊ�˰ѱ���ɫ�ŵ���ɫ�������ֲ���ʾ��colorbar���/20����������ɫ������20����ɫ
                    if (FPT_COLOR_TYPE=='C') %�������ɫͼ������Ҫ��colorbar
                        hc=colorbar;
                        set(hc,'ylim',[0 35],'FontName','times new roman','FontSize',12); %����colorbar����ʾ��Χ
                    end
                    
                    
                    eval(['[cc,hh]=contour(xr,yr,SPT_Comp_Data(:,:,' num2str(SPT_LAYER(k)) '),SPT_contourscale);']);
                    clabel(cc,hh,'FontSize',10,'FontName','Times new roman');
                    set(hh,'EdgeColor',[0.3 0.3 0.3],'linewidth',1) %���õ�ֵ�ߵ���ɫ�ʹ�ϸ
 
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
                    
                    %---------���ư��ߡ���------------------------
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
                    
                    if (FPT_COORDINATE=='BL') %����ͼ���������ѡȡ������Ӧ��label
                        xlabel('Longitude (\circE)');
                        ylabel('Latitude (\circN)');
                    elseif (FPT_COORDINATE=='XY')
                        xlabel('Distance (km)');
                        ylabel('Distance (km)');
                    end
                    
                    eval(['print(gcf,''-dpng'',''' fullfile(SPT_Pic_Path,SPT_Comp_Filename) '_' num2str(SPT_LAYER(k)) ''');']); %����ͼ��
                    
                end
            end
        end      
    end 
    %------------- End Salinity Field Distribution Drawing ----------------
    
    
    %----------- Residual Salinity Field Distribution Drawing -------------
    if (LOG_RSPT=='T') %�����Ҫ���ζȵ�ƽ��ͼ
        disp('Residual Salinity Field Distribution Drawing ...')
        RSPT_Pic_Path=fullfile(OUT_DIRE,'residual_distri\field\salinity'); %���ζ�ƽ��ͼ���·��
        [s,mess,messid]=mkdir(RSPT_Pic_Path); %�����ζ�ƽ��ͼ�Ĵ���ļ���
        RSPT_Comp_Filepath=fullfile(IN_DIRE,'resi_flux'); %ģʽ�����ζȳ��������ļ�·��
        
        RSPT_LAYER=str2num(RSPT_LAYER); %��������ζȻ�ͼ�������ַ�����ת��Ϊ��������
        Num_RSPT_LAYER=length(RSPT_LAYER); %�ܹ���Ҫ��������ε�ͼ��
        RSPT_contourscale=[-0.5 0 0.5 1 3 5 10 15 20 25 30 35]; %��Ҫ���ĵ�ֵ��
        
        RSPT_Comp_Fileinfo=dir(fullfile(RSPT_Comp_Filepath,'resi_salt_flux_3d_*.out')); %ģʽ�������ļ���Ϣ�������ļ������޸�ʱ�䣬��С���Ƿ�Ϊ�ļ�·����
        
        
        for i=1:length(RSPT_Comp_Fileinfo) %������Ҫ��ͼ����������ѭ��
            
            %--------------------- Residual Salinity Field Distribution Data Reading -----------------------
            RSPT_Comp_Filename=RSPT_Comp_Fileinfo(i).name; %ģʽ������������ζ��ļ�           
            eval(['fid_RSPT_Comp = fopen(''' fullfile(RSPT_Comp_Filepath,RSPT_Comp_Filename) ''',''r'');']); %�����ζ��ļ�
            RSPT_Comp_Data=nan*ones(IM,JM,KB);
            Note=fscanf(fid_RSPT_Comp,'%s',7);
            RSPT_B=fscanf(fid_RSPT_Comp,'%f',1);
            Note=fscanf(fid_RSPT_Comp,'%s',2);
            RSPT_E=fscanf(fid_RSPT_Comp,'%f',1);
            Note=fgetl(fid_RSPT_Comp);
            Note=fgetl(fid_RSPT_Comp);
            disp(['  Residual Salinty Distribution from ' num2str(RSPT_B) ' to ' num2str(RSPT_E) ' Hour Drawing'])
            
            while(~feof(fid_RSPT_Comp)) %�ж��Ƿ�����ļ�β��
                I=fscanf(fid_RSPT_Comp,'%d',1); %��I
                J=fscanf(fid_RSPT_Comp,'%d',1); %��J
                Note=fscanf(fid_RSPT_Comp,'%f',2); %��X,Y
                RSPT_Comp_Data(I,J,1:KB-1)=fscanf(fid_RSPT_Comp,'%f',5); %��X,Yλ��kbm1����ζ�ֵ
                Note=fgetl(fid_RSPT_Comp); %����ĩβ�Ļس�
            end
            
            fclose(fid_RSPT_Comp); %�رռ�����������ζ��ļ�
            %------------------- End Residual Salinity Field Distribution Data Reading -----------------------
            
            for i1=1:IM
                for j1=1:JM
                    if (h(i1,j1)<RS_MIN_DEP) %С��RS_MIN_DEPλ�ò���
                        RSPT_Comp_Data(i1,j1,:)=NaN; 
                        
                        for k=1:KB-1
                            if RSPT_Comp_Data(i1,j1,k)<0
                                RSPT_Comp_Data(i1,j1,k)=0;
                            end
                        end
                    end
                end
            end  
            
            
            
            
            Pic_time_b=datevec(Init_day+RSPT_B/24); %����ģʽ���õ���ʼʱ��������ζȳ�ͳ�Ƶ���ʼʱ��
            Pic_time_e=datevec(Init_day+RSPT_E/24); %����ģʽ���õ���ʼʱ��������ζȳ�ͳ�ƵĽ���ʱ��
            
            
            for k=1:Num_RSPT_LAYER
                close(figure(1));
                disp(['    Layer ' num2str(RSPT_LAYER(k))])
                if (FPT_COLOR_TYPE=='C') %���ݻ�ͼ����ɫ���ͣ�ѡ����ɫ���ļ�
                    map=load ('Colormap\cm_s_C.dat');
                elseif (FPT_COLOR_TYPE=='G')
                    map=load ('Colormap\cm_s_G.dat');
                end
                map=[bgcolor;map]; %��ӱ�����ɫ�����ζ�Ϊnan��λ�ð�����ɫ��
                colormap(map);  %������ɫ��
                
                
                set(gca,'box','on','FontName','times new roman','FontSize',12);
                set(gca,'color',bgcolor);
                set(gca,'xlim',[FPT_XMIN FPT_XMAX],'xtick',linspace(FPT_XMIN,FPT_XMAX,4));
                set(gca,'ylim',[FPT_YMIN FPT_YMAX],'ytick',linspace(FPT_YMIN,FPT_YMAX,4));
                set(gcf,'inverthardcopy','off'); %����ͼƬʱ�������õ���ɫ�����Զ�����
                set(gcf,'color',[1 1 1]);
                hold on 
                
                
                warning off; %��contourfʱ���кܶྯ����Ϣ������Ӱ�컭ͼ�����Թص�warning
                eval(['contourf(xr,yr,RSPT_Comp_Data(:,:,' num2str(RSPT_LAYER(k)) '),RSPT_contourscale);']);
                shading flat; 
                caxis([-0.5 35]); %Ϊ�˰ѱ���ɫ�ŵ���ɫ�������ֲ���ʾ��colorbar���/20����������ɫ������20����ɫ
                if (FPT_COLOR_TYPE=='C') %�������ɫͼ������Ҫ��colorbar
                    hc=colorbar;
                    set(hc,'ylim',[0 35],'FontName','times new roman','FontSize',12); %����colorbar����ʾ��Χ
                end
                
                
                eval(['[cc,hh]=contour(xr,yr,RSPT_Comp_Data(:,:,' num2str(RSPT_LAYER(k)) '),RSPT_contourscale);']);
                clabel(cc,hh,'FontSize',10,'FontName','Times new roman');
                set(hh,'EdgeColor',[0.3 0.3 0.3],'linewidth',1) %���õ�ֵ�ߵ���ɫ�ʹ�ϸ
                
                
                
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
                
                %---------���ư��ߡ���------------------------
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
                
                if (FPT_COORDINATE=='BL') %����ͼ���������ѡȡ������Ӧ��label
                    xlabel('Longitude (\circE)');
                    ylabel('Latitude (\circN)');
                elseif (FPT_COORDINATE=='XY')
                    xlabel('Distance (km)');
                    ylabel('Distance (km)');
                end
                
                eval(['print(gcf,''-dpng'',''' fullfile(RSPT_Pic_Path,RSPT_Comp_Filename(1:end-4)) '_' num2str(RSPT_LAYER(k)) ''');']); %����ͼ��
                
            end
        end
    end 
    %-------- End Residual Salinity Field Distribution Drawing -------------    
    
    
    %------------- Sediment Field Distribution Drawing -------------------
    if (LOG_SED=='T') %�����Ҫ����ɳ����ƽ��ͼ
        disp('Sediment Field Distribution Drawing ...')
        SED_Pic_Path=fullfile(OUT_DIRE,'field_distri\sediment'); %��ɳ��ƽ��ͼ���·��
        [s,mess,messid]=mkdir(SED_Pic_Path); %������ɳ��ƽ��ͼ�Ĵ���ļ���
        EPT_Comp_Filepath=fullfile(IN_DIRE,'field_distri\elevation'); %ģʽ��ˮλ���������ļ�·��
        SED_Comp_Filepath=fullfile(IN_DIRE,'field_distri\sediment'); %ģʽ�ĺ�ɳ�����������ļ�·��
        
        
        SED_TIME=str2num(SED_TIME); %������ĺ�ɳ����ͼʱ�̵��ַ�����ת��Ϊ��������
        Num_SED_TIME=length(SED_TIME); %�ܹ���Ҫ������ʱ�̵�ͼ��
        NSED_LAYER=str2num(SED_LAYER); %������ĺ�ɳ����ͼ�������ַ�����ת��Ϊ��������
        Num_SED_LAYER=length(NSED_LAYER); %�ܹ���Ҫ��������ε�ͼ��
        SED_contourscale=[-100 0 0.1 0.5 1 1.5 2 3 5]; %��Ҫ���ĵ�ֵ��
            
        
        
        for i=1:Num_SED_TIME %������Ҫ��ͼ����������ѭ��
            
            disp(['  The ' num2str(SED_TIME(i)) ' Hours Drawing'])
            
            %------ Elevation Field Distribution Data Reading ,For Check Tideflat Grid ----------
            EPT_Comp_Filename_SN=SED_TIME(i)/(N_FPT/3600); %����ģʽ�����ʱ������������Ҫ��ͼʱ�̶�Ӧ�����
            EPT_Comp_Filename=sprintf('el_field_%06.6d',EPT_Comp_Filename_SN); %ģʽ���������ˮλ�ļ�
            SED_Comp_Filename_SN=SED_TIME(i)/(N_FPT/3600); %����ģʽ�����ʱ������������Ҫ��ͼʱ�̶�Ӧ�����
            SED_Comp_Filename=sprintf('sed_field_%06.6d',SED_Comp_Filename_SN); %ģʽ���������ˮλ�ļ�           
            
            if ~(exist(fullfile(EPT_Comp_Filepath,EPT_Comp_Filename))&&exist(fullfile(SED_Comp_Filepath,SED_Comp_Filename))) %�ж�ˮλ���ļ��ͺ�ɳ�����ļ��Ƿ񶼴��ڣ��������������ʾ��Ϣ
                disp(['    ' EPT_Comp_Filename ' or ' SED_Comp_Filename ' dose not exist, drawing failed, please check the setting!']) %�����ļ������ڣ������ʾ��Ϣ
            else %�����ļ����ڣ����ͼ
                eval(['fid_EPT_Comp = fopen(''' fullfile(EPT_Comp_Filepath,EPT_Comp_Filename) ''',''r'',''b'');']); %���ļ������ڸö������ļ����á�BIG_ENDIAN����װ�����Դ���ѡ�������b��
                brecord=fread(fid_EPT_Comp,1,'integer*4'); %����¼��С��Ϣ
                thour=fread(fid_EPT_Comp,1,'real*4');
                EPT_Comp_Data=fread(fid_EPT_Comp,IM*JM,'real*4');
                erecord=fread(fid_EPT_Comp,1,'integer*4'); %����¼��С��Ϣ
                
                if (brecord~=erecord) %���һ����¼��ǰ������¼�Ĵ�С��һ�£����˳�����
                    disp(['Error in reading file ' fullfile(EPT_Comp_Filepath,EPT_Comp_Filename) ', please check!'])
                    break
                end
                EPT_Comp_Data=reshape(EPT_Comp_Data,IM,JM); %��ˮλ����ת��Ϊ������ʽ
                fclose(fid_EPT_Comp); %�رռ��������ˮλ�ļ�
                %---------------- End Elevation Field Distribution Data Reading ----------------------
                
                
                
                %--------------------- Sediment Field Distribution Data Reading -----------------------
                eval(['fid_SED_Comp = fopen(''' fullfile(SED_Comp_Filepath,SED_Comp_Filename) ''',''r'',''b'');']); %���ļ������ڸö������ļ����á�BIG_ENDIAN����װ�����Դ���ѡ�������b��
                brecord=fread(fid_SED_Comp,1,'integer*4'); %����¼��С��Ϣ
                thour=fread(fid_SED_Comp,1,'real*4');
                SED_Comp_Data=fread(fid_SED_Comp,IM*JM*KB,'real*4');
                erecord=fread(fid_SED_Comp,1,'integer*4'); %����¼��С��Ϣ
                
                if (brecord~=erecord) %���һ����¼��ǰ������¼�Ĵ�С��һ�£����˳�����
                    disp(['Error in reading file ' fullfile(SED_Comp_Filepath,SED_Comp_Filename) ', please check!'])
                    break
                end
                SED_Comp_Data=reshape(SED_Comp_Data,IM,JM,KB); %����ɳ��ת��Ϊ������ʽ
                fclose(fid_SED_Comp); %�رռ��������ˮλ�ļ�
                %------------------- End Sediment Field Distribution Data Reading -----------------------
                
                
                for i1=1:IM
                    for j1=1:JM
                        if (h(i1,j1)+EPT_Comp_Data(i1,j1)<DMIN||h(i1,j1)<-10) %����̲��½�ش��ĺ�ɳ����ֵΪnan
                            SED_Comp_Data(i1,j1,:)=NaN;
                        end
                        
                        for k=1:KB-1
                            if SED_Comp_Data(i1,j1,k)<0
                                SED_Comp_Data(i1,j1,k)=0;
                            end
                        end
                    end
                end  
                
                
                
                Pic_time=datevec(Init_day+thour/24); %����ģʽ���õ���ʼʱ��͸����ݼ�¼��thour����ʵ�ʵ�ʱ�䣬Pic_time���� �ꡢ�¡��ա�ʱ���֡��� ��Ϣ
                
                
                for k=1:Num_SED_LAYER
                    close(figure(1));
                    disp(['    Layer ' num2str(SED_LAYER(k))])
                    if (FPT_COLOR_TYPE=='C') %���ݻ�ͼ����ɫ���ͣ�ѡ����ɫ���ļ�
                        map=load ('Colormap\cm_sed2_C.dat');
                    elseif (FPT_COLOR_TYPE=='G')
                        map=load ('Colormap\cm_sed_G.dat');
                    end
                     map=[bgcolor;map]; %��ӱ�����ɫ���ú�ɳΪnan��λ�ð�����ɫ��
                    colormap(map);  %������ɫ��
                    
                    
                    set(gca,'box','on','FontName','times new roman','FontSize',12);
                    set(gca,'color',bgcolor);
                    set(gca,'xlim',[FPT_XMIN FPT_XMAX],'xtick',linspace(FPT_XMIN,FPT_XMAX,4));
                    set(gca,'ylim',[FPT_YMIN FPT_YMAX],'ytick',linspace(FPT_YMIN,FPT_YMAX,4));
                    set(gcf,'inverthardcopy','off'); %����ͼƬʱ�������õ���ɫ�����Զ�����
                    set(gcf,'color',[1 1 1]);
                    hold on                    
                    
                    warning off; %��contourfʱ���кܶྯ����Ϣ������Ӱ�컭ͼ�����Թص�warning
                    eval(['contourf(xr,yr,SED_Comp_Data(:,:,' num2str(SED_LAYER(k)) '),SED_contourscale);']);
                    shading flat; 
                                
                    cmin=0;
                    cmax=4;
                     caxis([cmin-(cmax-cmin)/(length(map)-1) cmax]); %Ϊ�˰ѱ���ɫ�ŵ���ɫ�������ֲ���ʾ��colorbar���/20����������ɫ������20����ɫ
                    if (FPT_COLOR_TYPE=='C') %�������ɫͼ������Ҫ��colorbar
                        hc=colorbar;
                         set(hc,'ylim',[cmin cmax],'FontName','times new roman','FontSize',12); %����colorbar����ʾ��Χ
                    end
  
      
                    
                    eval(['[cc,hh]=contour(xr,yr,SED_Comp_Data(:,:,' num2str(SED_LAYER(k)) '),SED_contourscale);']);
                    clabel(cc,hh,'FontSize',10,'FontName','Times new roman');
                    set(hh,'EdgeColor',[0.3 0.3 0.3],'linewidth',1) %���õ�ֵ�ߵ���ɫ�ʹ�ϸ
                    
                    
                    
                    % 					if (SED_LAYER(k)==1)
                    % 						SED_LAYER_NAME='��';
                    % 					elseif (SED_LAYER(k)==5)
                    % 						SED_LAYER_NAME='��';
                    % 					else
                    % 						SED_LAYER_NAME=['��' num2str(SED_LAYER(k))];
                    % 					end
                    % 					title([num2str(Pic_time(1)) '��' num2str(Pic_time(2)) '��' num2str(Pic_time(3)) '��' num2str(Pic_time(4)) 'ʱ ' SED_LAYER_NAME '�㺬ɳ���ֲ�'])
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
                    
                    %---------���ư��ߡ���------------------------
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
                    
                    if (FPT_COORDINATE=='BL') %����ͼ���������ѡȡ������Ӧ��label
                        xlabel('Longitude (\circE)');
                        ylabel('Latitude (\circN)');
                    elseif (FPT_COORDINATE=='XY')
                        xlabel('Distance (km)');
                        ylabel('Distance (km)');
                    end
                    
                    eval(['print(gcf,''-dpng'',''' fullfile(SED_Pic_Path,SED_Comp_Filename) '_' num2str(SED_LAYER(k)) ''');']); %����ͼ��
                    
                end
            end
        end      
    end 
    %------------- End Sediment Field Distribution Drawing ----------------   
    
    %----------- Residual Sedimnet Field Distribution Drawing -------------
      
    if (LOG_RSED=='T') %�����Ҫ����ɳ����ƽ��ͼ
        
        disp('Residual sediment Field Distribution Drawing ...')
        RSED_Pic_Path=fullfile(OUT_DIRE,'residual_distri\field\sediment'); %�ຬɳ��ƽ��ͼ���·��
        [s,mess,messid]=mkdir(RSED_Pic_Path); %������ɳ��ƽ��ͼ�Ĵ���ļ���
        RSED_Comp_Filepath=fullfile(IN_DIRE,'resi_flux'); %ģʽ���ຬɳ�����������ļ�·��
        
 
        RSED_contourscale=[-0.1 0 0.1 0.5 1 1.5 2 5]; %��Ҫ���ĵ�ֵ��
        
        RSED_Comp_Fileinfo=dir(fullfile(RSED_Comp_Filepath,'resi_sed_flux_3d_*.out')); %ģʽ�������ļ���Ϣ�������ļ������޸�ʱ�䣬��С���Ƿ�Ϊ�ļ�·����
        
        
        for i=1:length(RSED_Comp_Fileinfo) %������Ҫ��ͼ����������ѭ��
            
            %--------------------- Residual sediment Field Distribution Data Reading -----------------------
            RSED_Comp_Filename=RSED_Comp_Fileinfo(i).name; %ģʽ����������ຬɳ���ļ�           
            eval(['fid_RSED_Comp = fopen(''' fullfile(RSED_Comp_Filepath,RSED_Comp_Filename) ''',''r'');']); %���ຬɳ���ļ�
            RSED_Comp_Data=nan*ones(IM,JM,KB-1);
            Note=fscanf(fid_RSED_Comp,'%s',7);
            RSED_B=fscanf(fid_RSED_Comp,'%f',1);
            Note=fscanf(fid_RSED_Comp,'%s',2);
            RSED_E=fscanf(fid_RSED_Comp,'%f',1);
            Note=fgetl(fid_RSED_Comp);
            Note=fgetl(fid_RSED_Comp);
            disp(['  Residual Sediment Distribution from ' num2str(RSED_B) ' to ' num2str(RSED_E) ' Hour Drawing'])
            
            while(~feof(fid_RSED_Comp)) %�ж��Ƿ�����ļ�β��
                I=fscanf(fid_RSED_Comp,'%d',1); %��I
                J=fscanf(fid_RSED_Comp,'%d',1); %��J
                Note=fscanf(fid_RSED_Comp,'%f',2); %��X,Y
                RSED_Comp_Data(I,J,1:KB-1)=fscanf(fid_RSED_Comp,'%f',KB-1); %��X,Yλ��kbm1��ĺ�ɳ��ֵ
                Note=fgetl(fid_RSED_Comp); %����ĩβ�Ļس�
            end
            
            fclose(fid_RSED_Comp); %�رռ���������ຬɳ���ļ�
            %------------------- End Residual sediment Field Distribution Data Reading -----------------------
            
            for i1=1:IM
                for j1=1:JM
                    if (h(i1,j1)<RS_MIN_DEP) %С��RS_MIN_DEPλ�ò���
                        RSED_Comp_Data(i1,j1,:)=NaN; 
                        
                        for k=1:KB-1
                            if RSED_Comp_Data(i1,j1,k)<0
                                RSED_Comp_Data(i1,j1,k)=0;
                            end
                        end
                    end
                end
            end  
            
            
            
            
            Pic_time_b=datevec(Init_day+RSED_B/24); %����ģʽ���õ���ʼʱ������ຬɳ����ͳ�Ƶ���ʼʱ��
            Pic_time_e=datevec(Init_day+RSED_E/24); %����ģʽ���õ���ʼʱ������ຬɳ����ͳ�ƵĽ���ʱ��
           
        RSED_LAYER=str2num(RSED_LAYER); %������ĺ�ɳ����ͼ�������ַ�����ת��Ϊ��������
        Num_RSED_LAYER=length(RSED_LAYER); %�ܹ���Ҫ��������ε�ͼ��
            
            for k=1:Num_RSED_LAYER
                close(figure(1));
                disp(['    Layer ' num2str(RSED_LAYER(k))])
                if (FPT_COLOR_TYPE=='C') %���ݻ�ͼ����ɫ���ͣ�ѡ����ɫ���ļ�
                    map=load ('Colormap\cm_sed2_C.dat');
                elseif (FPT_COLOR_TYPE=='G')
                    map=load ('Colormap\cm_sed_G.dat');
                end
                map=[bgcolor;map]; %��ӱ�����ɫ���ú�ɳ��Ϊnan��λ�ð�����ɫ��
                colormap(map);  %������ɫ��
                
                set(gcf,'position',[50 50 1000 600]);
                x0=0.05;  y0=0.1; width=0.85; height=0.8;
                    
                 subplot('position',[x0 y0 width height]);
                 axis equal;
                set(gca,'box','on','FontName','times new roman','FontSize',12);
                set(gca,'color',bgcolor);
                set(gca,'xlim',[FPT_XMIN FPT_XMAX],'xtick',linspace(FPT_XMIN,FPT_XMAX,4));
                set(gca,'ylim',[FPT_YMIN FPT_YMAX],'ytick',linspace(FPT_YMIN,FPT_YMAX,4));
                set(gcf,'inverthardcopy','off'); %����ͼƬʱ�������õ���ɫ�����Զ�����
                set(gcf,'color',[1 1 1]);

                hold on 
                
                
                warning off; %��contourfʱ���кܶྯ����Ϣ������Ӱ�컭ͼ�����Թص�warning
                eval(['contourf(xr,yr,RSED_Comp_Data(:,:,' num2str(RSED_LAYER(k)) '),RSED_contourscale);']);
                shading flat; 
%                 caxis([-0.1 5]); %Ϊ�˰ѱ���ɫ�ŵ���ɫ�������ֲ���ʾ��colorbar���/20����������ɫ������20����ɫ
%                 if (FPT_COLOR_TYPE=='C') %�������ɫͼ������Ҫ��colorbar
%                     hc=colorbar;
%                     set(hc,'ylim',[0 5],'FontName','times new roman','FontSize',12); %����colorbar����ʾ��Χ
%                 end
                cmin=0;
                cmax=5;
                caxis([cmin-(cmax-cmin)/(length(map)-1) cmax]); %Ϊ�˰ѱ���ɫ�ŵ���ɫ�������ֲ���ʾ��colorbar���/20����������ɫ������20����ɫ
                if (FPT_COLOR_TYPE=='C') %�������ɫͼ������Ҫ��colorbar
                    hc=colorbar;
                    set(hc,'ylim',[cmin cmax],'FontName','times new roman','FontSize',12); %����colorbar����ʾ��Χ
                end
                
                
                eval(['[cc,hh]=contour(xr,yr,RSED_Comp_Data(:,:,' num2str(RSED_LAYER(k)) '),RSED_contourscale);']);
                clabel(cc,hh,'FontSize',10,'FontName','Times new roman');
                set(hh,'EdgeColor',[0.3 0.3 0.3],'linewidth',1) %���õ�ֵ�ߵ���ɫ�ʹ�ϸ
                
                
                
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
                
                %---------���ư��ߡ���------------------------
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
                
                if (FPT_COORDINATE=='BL') %����ͼ���������ѡȡ������Ӧ��label
                    xlabel('Longitude (\circE)');
                    ylabel('Latitude (\circN)');
                elseif (FPT_COORDINATE=='XY')
                    xlabel('Distance (km)');
                    ylabel('Distance (km)');
                end
                
                eval(['print(gcf,''-dpng'',''' fullfile(RSED_Pic_Path,RSED_Comp_Filename(1:end-4)) '_' num2str(RSED_LAYER(k)) ''');']); %����ͼ��
                
            end
        end
    end 
    %-------- End Residual sediment Field Distribution Drawing -------------       
    
    %------------- Bottom Shear Stress Field Distribution Drawing -------------------
    if (LOG_TAU=='T') %�����Ҫ����Ӧ����ƽ��ͼ
        disp('Bottom Shear Stress Field Distribution Drawing ...')
        TAU_Pic_Path=fullfile(OUT_DIRE,'field_distri\stress'); %��Ӧ��ƽ��ͼ���·��
        [s,mess,messid]=mkdir(TAU_Pic_Path); %������Ӧ��ƽ��ͼ�Ĵ���ļ���
        EPT_Comp_Filepath=fullfile(IN_DIRE,'field_distri\elevation'); %ģʽ��ˮλ���������ļ�·��
        TAU_Comp_Filepath=fullfile(IN_DIRE,'field_distri\stress'); %ģʽ����Ӧ�����������ļ�·��
        
        
        TAU_TIME=str2num(TAU_TIME); %���������Ӧ����ͼʱ�̵��ַ�����ת��Ϊ��������
        Num_TAU_TIME=length(TAU_TIME); %�ܹ���Ҫ������ʱ�̵�ͼ��
        %        TAU_LAYER=str2num(TAU_LAYER); %���������Ӧ����ͼ�������ַ�����ת��Ϊ��������
        %        Num_TAU_LAYER=length(TAU_LAYER); %�ܹ���Ҫ��������ε�ͼ��
        TAU_contourscale=[-0.1 0 0.5 1 3 5]; %��Ҫ���ĵ�ֵ��
        
        for i=1:Num_TAU_TIME %������Ҫ��ͼ����������ѭ��
            
            disp(['  The ' num2str(TAU_TIME(i)) ' Hours Drawing'])
            
            %------ Elevation Field Distribution Data Reading ,For Check Tideflat Grid ----------
            EPT_Comp_Filename_SN=TAU_TIME(i)/(N_FPT/3600); %����ģʽ�����ʱ������������Ҫ��ͼʱ�̶�Ӧ�����
            EPT_Comp_Filename=sprintf('el_field_%06.6d',EPT_Comp_Filename_SN); %ģʽ���������ˮλ�ļ�
            TAU_Comp_Filename_SN=TAU_TIME(i)/(N_FPT/3600); %����ģʽ�����ʱ������������Ҫ��ͼʱ�̶�Ӧ�����
            TAU_Comp_Filename=sprintf('tau_field_%06.6d',TAU_Comp_Filename_SN); %ģʽ���������ˮλ�ļ�           
            
            if ~(exist(fullfile(EPT_Comp_Filepath,EPT_Comp_Filename))&&exist(fullfile(TAU_Comp_Filepath,TAU_Comp_Filename))) %�ж�ˮλ���ļ����ζȳ��ļ��Ƿ񶼴��ڣ��������������ʾ��Ϣ
                disp(['    ' EPT_Comp_Filename ' or ' TAU_Comp_Filename ' dose not exist, drawing failed, please check the setting!']) %�����ļ������ڣ������ʾ��Ϣ
            else %�����ļ����ڣ����ͼ
                eval(['fid_EPT_Comp = fopen(''' fullfile(EPT_Comp_Filepath,EPT_Comp_Filename) ''',''r'',''b'');']); %���ļ������ڸö������ļ����á�BIG_ENDIAN����װ�����Դ���ѡ�������b��
                brecord=fread(fid_EPT_Comp,1,'integer*4'); %����¼��С��Ϣ
                thour=fread(fid_EPT_Comp,1,'real*4');
                EPT_Comp_Data=fread(fid_EPT_Comp,IM*JM,'real*4');
                erecord=fread(fid_EPT_Comp,1,'integer*4'); %����¼��С��Ϣ
                
                if (brecord~=erecord) %���һ����¼��ǰ������¼�Ĵ�С��һ�£����˳�����
                    disp(['Error in reading file ' fullfile(EPT_Comp_Filepath,EPT_Comp_Filename) ', please check!'])
                    break
                end
                EPT_Comp_Data=reshape(EPT_Comp_Data,IM,JM); %��ˮλ����ת��Ϊ������ʽ
                fclose(fid_EPT_Comp); %�رռ��������ˮλ�ļ�
                %---------------- End Elevation Field Distribution Data Reading ----------------------
                
                
                
                %--------------------- TAUiment Field Distribution Data Reading -----------------------
                eval(['fid_TAU_Comp = fopen(''' fullfile(TAU_Comp_Filepath,TAU_Comp_Filename) ''',''r'',''b'');']); %���ļ������ڸö������ļ����á�BIG_ENDIAN����װ�����Դ���ѡ�������b��
                brecord=fread(fid_TAU_Comp,1,'integer*4'); %����¼��С��Ϣ
                thour=fread(fid_TAU_Comp,1,'real*4');
                TAU_Comp_Data=fread(fid_TAU_Comp,IM*JM,'real*4');
                erecord=fread(fid_TAU_Comp,1,'integer*4'); %����¼��С��Ϣ
                
                if (brecord~=erecord) %���һ����¼��ǰ������¼�Ĵ�С��һ�£����˳�����
                    disp(['Error in reading file ' fullfile(TAU_Comp_Filepath,TAU_Comp_Filename) ', please check!'])
                    break
                end
                TAU_Comp_Data=reshape(TAU_Comp_Data,IM,JM); %����Ӧ��ת��Ϊ������ʽ
                fclose(fid_TAU_Comp); %�رռ��������ˮλ�ļ�
                %------------------- End Salinity Field Distribution Data Reading -----------------------
                
                
                for i1=1:IM
                    for j1=1:JM
                        if (h(i1,j1)+EPT_Comp_Data(i1,j1)<DMIN||h(i1,j1)<-10) %����̲��½�ش�����Ӧ����ֵΪnan
                            TAU_Comp_Data(i1,j1,:)=NaN;
                        end
                    end
                end  
                
                
                
                Pic_time=datevec(Init_day+thour/24); %����ģʽ���õ���ʼʱ��͸����ݼ�¼��thour����ʵ�ʵ�ʱ�䣬Pic_time���� �ꡢ�¡��ա�ʱ���֡��� ��Ϣ
                
                
                %				for k=1:Num_TAU_LAYER
                close(figure(1));
                %					disp(['    Layer ' num2str(TAU_LAYER(k))])
                if (FPT_COLOR_TYPE=='C') %���ݻ�ͼ����ɫ���ͣ�ѡ����ɫ���ļ�
                    map=load ('Colormap\cm_sed_C.dat');
                elseif (FPT_COLOR_TYPE=='G')
                    map=load ('Colormap\cm_sed_G.dat');
                end
                map=[bgcolor;map]; %��ӱ�����ɫ�����ζ�Ϊnan��λ�ð�����ɫ��
                colormap(map);  %������ɫ��
                
                
                set(gca,'box','on','FontName','times new roman','FontSize',12);
                set(gca,'color',bgcolor);
                set(gca,'xlim',[FPT_XMIN FPT_XMAX],'xtick',linspace(FPT_XMIN,FPT_XMAX,4));
                set(gca,'ylim',[FPT_YMIN FPT_YMAX],'ytick',linspace(FPT_YMIN,FPT_YMAX,4));
                set(gcf,'inverthardcopy','off'); %����ͼƬʱ�������õ���ɫ�����Զ�����
                set(gcf,'color',[1 1 1]);
                hold on 
                
                
                warning off; %��contourfʱ���кܶྯ����Ϣ������Ӱ�컭ͼ�����Թص�warning
                eval(['contourf(xr,yr,TAU_Comp_Data(:,:),TAU_contourscale);']);
                shading flat; 
                caxis([-0.1 5]); %Ϊ�˰ѱ���ɫ�ŵ���ɫ�������ֲ���ʾ��colorbar���/20����������ɫ������20����ɫ
                if (FPT_COLOR_TYPE=='C') %�������ɫͼ������Ҫ��colorbar
                    hc=colorbar;
                    set(hc,'ylim',[0 5],'FontName','times new roman','FontSize',12,'YColor',[0 0 0]); %����colorbar����ʾ��Χ
                end
                
                
                eval(['[cc,hh]=contour(xr,yr,TAU_Comp_Data(:,:),TAU_contourscale);']);
                clabel(cc,hh,'FontSize',10,'FontName','Times new roman');
                set(hh,'EdgeColor',[0.3 0.3 0.3],'linewidth',1) %���õ�ֵ�ߵ���ɫ�ʹ�ϸ
                

                hold on  
                
                %---------���ư��ߡ���------------------------
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
                
                if (FPT_COORDINATE=='BL') %����ͼ���������ѡȡ������Ӧ��label
                    xlabel('Longitude (\circE)');
                    ylabel('Latitude (\circN)');
                elseif (FPT_COORDINATE=='XY')
                    xlabel('Distance (km)');
                    ylabel('Distance (km)');
                end
                
                eval(['print(gcf,''-dpng'',''' fullfile(TAU_Pic_Path,TAU_Comp_Filename) ''');']); %����ͼ��
                
            end
        end      
    end
    %------------- End Bottom Shear Stress Field Distribution Drawing ----------------



