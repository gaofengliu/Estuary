function re=plain_interp(sec_x,sec_y,var,im,jm,kb,xr,yr,fsm,h,h1,h2,extraval) 
%ƽ��˫���Բ�ֵ���㷨��Ecomһ��
%ע��sec_x,sec_yΪ��Ҫ��ֵ��λ�ã�varΪ������ֵ�ı�����im��jm��kbΪ�����С��xr��yrΪ����ʵ��λ�ã�fsmΪ��ʪ�����־
%hΪˮ�h1��h2Ϊ���񳤶ȣ�extravalΪ����½��ʱ��ȱʡֵ
num_point=length(sec_x); %�����м�������Ҫ��ֵ
for n=1:num_point
    %===========������һ����������ı��=============
    %ע��½�������У��󲿷�xr��yrΪ0
    x0=sec_x(n);
    y0=sec_y(n);
    dpmin=10000000;
    for i=1:im %�����ֵ���������������
        for j=1:jm
            % 			if h(i,j)>=-10
            dm=sqrt((x0-xr(i,j))^2+(y0-yr(i,j))^2);
            if dm<=dpmin
                dpmin=dm;
                i0=i;
                j0=j;
            end
            % 			end
        end
    end
    % 	sec_i(n)=i0;
    % 	sec_j(n)=j0;
    %======================End=====================
    
    %===========˫���Բ�ֵ=============
    r1=x0-xr(i0,j0);
    r2=y0-yr(i0,j0);
    d1=xr(i0+1,j0)-xr(i0,j0);
    d2=yr(i0+1,j0)-yr(i0,j0);
    c1=xr(i0,j0+1)-xr(i0,j0);
    c2=yr(i0,j0+1)-yr(i0,j0);
    dd1=(r1*d1+r2*d2)/sqrt(d1^2+d2^2);
    cc1=(r1*c1+r2*c2)/sqrt(c1^2+c2^2);
    
    if dd1>=0 
        i1=i0;
    else
        i1=i0-1;
    end
    
    if cc1>=0
        j1=j0;
    else
        j1=j0-1;
    end
    
    
    
    r11=x0-xr(i1,j1);
    r12=y0-yr(i1,j1);
    
    d11=xr(i1+1,j1)-xr(i1,j1);
    d12=yr(i1+1,j1)-yr(i1,j1);
    
    rr1=((r11*d11+r12*d12)/sqrt(d11^2+d12^2));
    cx1=rr1/h1(i1+1,j1);
    
    if cx1<=0
        cx1=0;
    elseif cx1>=1
        cx1=1;
    end
    cx2=1-cx1;
    
    d11=xr(i1,j1+1)-xr(i1,j1);
    d12=yr(i1,j1+1)-yr(i1,j1);
    
    rr1=((r11*d11+r12*d12)/sqrt(d11^2+d12^2));
    
    cy1=rr1/h2(i1,j1+1);
    if cy1<=0
        cy1=0;
    elseif cy1>=1
        cy1=1;
    end
    cy2=1-cy1;
    
    if length(size(var))==3 %��ֵ����var�������ά���������ζȻ�����
        for k=1:kb-1
            for j=j1:j1+1
                for i=i1:i1+1
                    YT(i,j)=var(i,j,k);
                    if fsm(i,j)==0
                        YT(i,j)=extraval; %��½�ش���ֵΪȱʡֵ
                    end
                end
            end
            
            hh1=cx2*YT(i1,j1)+cx1*YT(i1+1,j1);
            
            if fsm(i1,j1)==0
                hh1=YT(i1+1,j1);
            end
            
            if fsm(i1+1,j1)==0
                hh1=YT(i1,j1);
            end
            
            hh2=cx2*YT(i1,j1+1)+cx1*YT(i1+1,j1+1);
            
            if fsm(i1,j1+1)==0
                hh2=YT(i1+1,j1+1);
            end
            
            if fsm(i1+1,j1+1)==0
                hh2=YT(i1,j1+1);
            end
            re(k,n)=cy2*hh1+cy1*hh2;
            
            if (fsm(i1,j1)+fsm(i1+1,j1))==0
                re(k,n)=hh2;
            end
            
            if (fsm(i1,j1+1)+fsm(i1+1,j1+1))==0
                re(k,n)=hh1;
            end
            
            if fsm(i0,j0)==0 %����������һ�����������Ϊ�ɣ���ô�����ϸõ��ֵ��Ϊnan
                re(k,n)=extraval;
            end   
        end
        
    elseif length(size(var))==2 %��ֵ����var����Ƕ�ά��������ˮ���ˮλ
        for j=j1:j1+1
            for i=i1:i1+1
                YT(i,j)=var(i,j);
            end
        end
        
        hh1=cx2*YT(i1,j1)+cx1*YT(i1+1,j1);
        
        if fsm(i1,j1)==0
            hh1=YT(i1+1,j1);
        end
        
        if fsm(i1+1,j1)==0
            hh1=YT(i1,j1);
        end
        
        hh2=cx2*YT(i1,j1+1)+cx1*YT(i1+1,j1+1);
        
        if fsm(i1,j1+1)==0
            hh2=YT(i1+1,j1+1);
        end
        
        if fsm(i1+1,j1+1)==0
            hh2=YT(i1,j1+1);
        end
        re(n)=cy2*hh1+cy1*hh2;
        
        if (fsm(i1,j1)+fsm(i1+1,j1))==0
            re(n)=hh2;
        end
        if (fsm(i1,j1+1)+fsm(i1+1,j1+1))==0
            re(n)=hh1;
        end       
        
        if h(i0,j0)<=-10 %����������һ��������������ˮ��С��-10��Ϊ½�أ�����ô�����ϸõ�Ϊ��Ϊ��½�أ�ˮ�ֵΪSEC_MAX_HIGHT��Ϊ�˻�ͼ��Ҫˮλֵ��ΪSEC_MAX_HIGHT
            re(n)=extraval;
        end   
    end
end