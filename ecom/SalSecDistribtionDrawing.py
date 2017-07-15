
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
%------------- Salinity Section Distribution Drawing -------------------
    if (LOG_SEC_S=='T') %�����Ҫ���ζȵĶ���ͼ
        disp('Salinity Section Distribution Drawing ...')
        S_Sec_Pic_Path=fullfile(OUT_DIRE,'section_distri\salinity'); %�ζȶ���ͼ���·����ָ����
        [s,mess,messid]=mkdir(S_Sec_Pic_Path); %�����ζȶ���ͼ�Ĵ���ļ��У�ָ����
        EPT_Comp_Filepath=fullfile(IN_DIRE,'field_distri\elevation'); %ģʽ��ˮλ���������ļ�·��
        SPT_Comp_Filepath=fullfile(IN_DIRE,'field_distri\salinity'); %ģʽ���ζȳ��������ļ�·��
        
        SEC_S_TIME=str2num(SEC_S_TIME); %��������ζȶ����ͼʱ�̵��ַ�����ת��Ϊ��������
        Num_SEC_S_TIME=length(SEC_S_TIME); %�ܹ���Ҫ������ʱ�̵�ͼ��
        SPT_contourscale=[-0.5 0 0.5 1 3 5 10 15 20 25 30 35]; %��Ҫ���ĵ�ֵ��
        
        %��������߳�
        for i=2:IM
            for j=2:JM
                h1(i,j)=sqrt((xr(i-1,j)-xr(i,j))^2+(yr(i-1,j)-yr(i,j))^2);
                h2(i,j)=sqrt((xr(i,j-1)-xr(i,j))^2+(yr(i,j-1)-yr(i,j))^2);
            end
        end
       
        if strcmp(SEC_P_TYPE,'SPECIFY') %������������ļ���ָ�������λ��
            for nsec=1:SEC_NUM
                eval(['SEC_CONTROL_POINTS_' num2str(nsec) '=str2num(SEC_CONTROL_POINTS_' num2str(nsec) ');']); %������Ķ���λ�ô��ַ�������Ϊ�������
            end
        end
        
        for i=1:Num_SEC_S_TIME %������Ҫ��ͼ����������ѭ��
            
            disp(['  The ' num2str(SEC_S_TIME(i)) ' Hours Drawing'])
            
            %------ Elevation Field Distribution Data Reading ,For Check Tideflat Grid ----------
            EPT_Comp_Filename_SN=SEC_S_TIME(i)/(N_FPT/3600); %����ģʽ�����ʱ������������Ҫ��ͼʱ�̶�Ӧ�����
            EPT_Comp_Filename=sprintf('el_field_%06.6d',EPT_Comp_Filename_SN); %ģʽ���������ˮλ�ļ�
            SPT_Comp_Filename_SN=SEC_S_TIME(i)/(N_FPT/3600); %����ģʽ�����ʱ������������Ҫ��ͼʱ�̶�Ӧ�����
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
                
                fsm=ones(IM,JM);
                for i1=1:IM
                    for j1=1:JM
                        if (h(i1,j1)+EPT_Comp_Data(i1,j1)<DMIN||h(i1,j1)<-10) %����̲��½�ش����ζȸ�ֵΪnan�͸�ʪ������
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
                
                
                
                Pic_time=datevec(Init_day+thour/24); %����ģʽ���õ���ʼʱ��͸����ݼ�¼��thour����ʵ�ʵ�ʱ�䣬Pic_time���� �ꡢ�¡��ա�ʱ���֡��� ��Ϣ
                
                
                close(figure(1));
                close(figure(2));
                figure(1),
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
                contourf(xr,yr,SPT_Comp_Data(:,:,1),SPT_contourscale);
                shading flat; 
                caxis([-0.5 35]); %Ϊ�˰ѱ���ɫ�ŵ���ɫ�������ֲ���ʾ��colorbar���/20����������ɫ������20����ɫ
                if (FPT_COLOR_TYPE=='C') %�������ɫͼ������Ҫ��colorbar
                    hc=colorbar;
                    set(hc,'ylim',[0 35],'FontName','times new roman','FontSize',12); %����colorbar����ʾ��Χ
                end
                
                [cc,hh]=contour(xr,yr,SPT_Comp_Data(:,:,1),SPT_contourscale);
                clabel(cc,hh,'FontSize',10,'FontName','Times new roman');
                set(hh,'EdgeColor',[0.3 0.3 0.3],'linewidth',1) %���õ�ֵ�ߵ���ɫ�ʹ�ϸ
                
                %  				title([num2str(Pic_time(1)) '��' num2str(Pic_time(2)) '��' num2str(Pic_time(3)) '��' num2str(Pic_time(4)) 'ʱ ����ζȷֲ�'])
                title([num2str(Pic_time(4)) ':00 ' mon(Pic_time(2),:) ' ' num2str(Pic_time(3)) ' ' num2str(Pic_time(1)) ' (GMT +8)   Surface'])
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
             
                %=================== ���ƶ���ͼ ======================
                if strcmp(SEC_P_TYPE,'SPECIFY') %�������λ��Ԥ��ָ�� 
                    
                    for nsec=1:SEC_NUM %���ն�����Ŀ����ѭ��
                        eval(['SEC_CONTROL_POINTS=SEC_CONTROL_POINTS_' num2str(nsec) ';']) %��������Ƶ�ͳһ��ֵ��SEC_CONTROL_POINTS��
                        figure(1),plot(SEC_CONTROL_POINTS(:,1),SEC_CONTROL_POINTS(:,2),'linewidth',1.5,'color','k'); %��ͼ��1�У��ζȱ��ƽ��ͼ���л��������λ��
                        hold on
                        text(SEC_CONTROL_POINTS(1,1),SEC_CONTROL_POINTS(1,2),num2str(nsec),'FontName','times new roman','FontSize',15);
                        hold on
                        
                        SEC_S_POSSITION_X=ones(0); %��ʼ�����棨���ݷֱ��ʽ����Ƶ����䣩��������������ֵΪ��
                        SEC_S_POSSITION_Y=ones(0); %��ʼ�����棨���ݷֱ��ʽ����Ƶ����䣩��������������ֵΪ��
                        for ii=1:length(SEC_CONTROL_POINTS)-1 %�Զ��������������Ƶ���д���
                            xs=SEC_CONTROL_POINTS(ii,1); %��������
                            xe=SEC_CONTROL_POINTS(ii+1,1); %�յ������
                            ys=SEC_CONTROL_POINTS(ii,2); %���������
                            ye=SEC_CONTROL_POINTS(ii+1,2); %�յ�������
                            dis=sqrt((xs-xe)^2+(ys-ye)^2); %�����յ�ľ���
                            num=ceil(dis/SEC_RESOLUTION); %����Ԥ���趨�ķֱ���ȷ��ÿ�������Ƶ��ּ���
                            dx=(xe-xs)/num; %x�᷽��ÿ�εĳ���
                            dy=(ye-ys)/num; %y�᷽��ÿ�εĳ���
                            
                            for jj=1:num
                                xi=xs+(jj-1)*dx; %���������ĺ�����
                                yi=ys+(jj-1)*dy; %����������������
                                SEC_S_POSSITION_X=[SEC_S_POSSITION_X xi]; %��ӵ�����������
                                SEC_S_POSSITION_Y=[SEC_S_POSSITION_Y yi]; %��ӵ�����������
                            end
                            
                        end
                        SEC_S_POSSITION_X=[SEC_S_POSSITION_X SEC_CONTROL_POINTS(end,1)]; %�����һ�����Ƶ���ӵ�����������
                        SEC_S_POSSITION_Y=[SEC_S_POSSITION_Y SEC_CONTROL_POINTS(end,2)]; %�����һ�����Ƶ���ӵ�����������
                        
                        sec_h=plain_interp(SEC_S_POSSITION_X,SEC_S_POSSITION_Y,h,IM,JM,KB,xr,yr,fsm,h,h1,h2,SEC_MAX_HIGHT); %��ˮ�����˫���Բ�ֵ��½�ش���Ԥ���趨�����ͼ����߸̶߳��壨SEC_MAX_HIGHT��
                        sec_el=plain_interp(SEC_S_POSSITION_X,SEC_S_POSSITION_Y,EPT_Comp_Data,IM,JM,KB,xr,yr,fsm,h,h1,h2,-SEC_MAX_HIGHT); %��ˮλ����˫���Բ�ֵ��½�ش���Ԥ���趨�����ͼ����߸̵߳��෴�����壨-SEC_MAX_HIGHT�� ע�����ǵ�ˮλ��ˮ����������෴
                        sec_s=plain_interp(SEC_S_POSSITION_X,SEC_S_POSSITION_Y,SPT_Comp_Data,IM,JM,KB,xr,yr,fsm,h,h1,h2,nan); %���ζ���˫���Բ�ֵ��½�ش���nan��ʾ
                        
                        
                        %ע�������ζȶ�������������룬Ϊ�˵õ�������ͼ�񣬴��������¶���ֵ���������
                        for kk=1:length(SEC_S_POSSITION_X)
                            distance(kk)=sqrt((SEC_S_POSSITION_X(kk)-SEC_S_POSSITION_X(1))^2+(SEC_S_POSSITION_Y(kk)-SEC_S_POSSITION_Y(1))^2); %�������������������ľ���
                        end
                        x_sec_s=ones(KB+1,1)*distance; %��Ч���ζȿ���KB-1������������������2����ֵ�����Դ���Ϊ7������
                        
                        for ii=1:length(SEC_S_POSSITION_X)
                            sec_d=sec_h(ii)+sec_el(ii); %��ˮ�����ˮ�����ˮλ
                            for kk=1:KB-1    
                                y_sec_s(kk,ii)=-sec_el(ii)+sec_d/(KB-1)*(kk-1)+0.5*sec_d/(KB-1); %���м�KB-1����Ч������ڵ�λ�ã���ȣ�
                            end
                        end
                        y_sec_s=[-sec_el;y_sec_s;sec_h]; %���λ��ȡˮλ���෴�����ײ�λ��ȡ��Ӧ��ˮ��
                        
                        
                        sec_s=[sec_s(1,:);sec_s;sec_s(KB-1,:)]; %�ζȵ�ֵ���������
                        
                        
                        %�˴���������Ҫ�ĵ�ֵ�߷�Χ��Ӧ��֤��0.5��ֵ�ߣ���ˮ������δ������ȫ�����Ŀǰ����ֵ��ʱ����Ĭ�ϣ�       
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
                        
                        set(gca,'color',[0 1 0],'FontName','times new roman','FontSize',12); %����ɫ��Ϊ�����һ������ɫ          
                        set(gca,'xlim',[0 max(max(x_sec_s))]); %����������Ϊ0��������
                        set(gca,'ylim',[SEC_MAX_HIGHT ceil(max(max(y_sec_s)))]); %����������Ϊ��߸߳���ˮ�����ֵ������ȡ��
                        set(gca,'YDir','reverse') %��y�ᷴ��
                        set(gcf,'inverthardcopy','off'); %����ͼƬʱ�������õ���ɫ�����Զ�����
                        set(gcf,'color',[1 1 1]); 
                        hold on;
                        
                        %                          [cc,hh]=contour(x_sec_s,y_sec_s,sec_s,contourscale,'k');
                        
                        [cc,hh]=contour(x_sec_s,y_sec_s,sec_s,'k');
                        hold on;
                        clabel(cc,hh,'fontsize',10,'fontname','Times new roman','color','k','rotation',0);
                        hold on;
                        
                        fill([distance(1) distance distance(end)],[SEC_MAX_HIGHT -sec_el SEC_MAX_HIGHT],'w') %��ˮ���������Ϊ��ɫ
                        hold on                 
                        %                         plot(distance,-sec_el,'w')
                        hold on
                        
                        fill([distance(1) distance distance(end)],[ceil(max(max(y_sec_s))) sec_h ceil(max(max(y_sec_s)))],'g') %�����ΰ���½�����Ϊ��ɫ
                        hold on
                        plot(distance,sec_h,'g') %���ʱ�������������ɫ��һ�£��ʰѵ���������Ϊ��ɫ�����ڸ�
                        hold on
                        xlabel('Distance (km)');
                        ylabel('Depth (m)')
                        %                       title([num2str(Pic_time(1)) '��' num2str(Pic_time(2)) '��' num2str(Pic_time(3)) '��' num2str(Pic_time(4)) 'ʱ ����' num2str(nsec) '�ζȷֲ�'])
                        title([num2str(Pic_time(4)) ':00 ' mon(Pic_time(2),:) ' ' num2str(Pic_time(3)) ' ' num2str(Pic_time(1)) ' (GMT +8)   Sec' num2str(nsec)])                         
                        
                        figure(2),eval(['print(gcf,''-dpng'',''' fullfile(S_Sec_Pic_Path,SPT_Comp_Filename) '_specify_sec' num2str(nsec) ''');']); %�������ͼ��
                        
                        clear SEC_S_POSSITION_X SEC_S_POSSITION_Y x_sec_s y_sec_s distance sec_h sec_el sec_s %������ֱ���
                        
                    end
                    figure(1),eval(['print(gcf,''-dpng'',''' fullfile(S_Sec_Pic_Path,SPT_Comp_Filename) '_specify_sec_check'');']); %����ƽ��ͼ��
                    
             end
        end      
    end 
    %------------- End Salinity Section Distribution Drawing ----------------  
