function re=plain_interp(sec_x,sec_y,var,im,jm,kb,xr,yr,fsm,h,h1,h2,extraval) 
%平面双线性插值，算法与Ecom一致
%注：sec_x,sec_y为需要插值的位置；var为用来插值的变量；im，jm，kb为网格大小；xr，yr为网格实际位置；fsm为干湿网格标志
%h为水深；h1，h2为网格长度；extraval为碰到陆地时的缺省值
num_point=length(sec_x); %计算有几个点需要插值
for n=1:num_point
    %===========断面上一点所在网格的编号=============
    %注：陆地网格中，大部分xr，yr为0
    x0=sec_x(n);
    y0=sec_y(n);
    dpmin=10000000;
    for i=1:im %求离插值点最近的网格坐标
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
    
    %===========双线性插值=============
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
    
    if length(size(var))==3 %插值变量var如果是三维向量，即盐度或流速
        for k=1:kb-1
            for j=j1:j1+1
                for i=i1:i1+1
                    YT(i,j)=var(i,j,k);
                    if fsm(i,j)==0
                        YT(i,j)=extraval; %在陆地处赋值为缺省值
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
            
            if fsm(i0,j0)==0 %如果离断面上一点最近的网格为干，那么断面上该点的值赋为nan
                re(k,n)=extraval;
            end   
        end
        
    elseif length(size(var))==2 %插值变量var如果是二维向量，即水深或水位
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
        
        if h(i0,j0)<=-10 %如果离断面上一点最近的网格点总水深小于-10（为陆地），那么断面上该点为认为是陆地，水深赋值为SEC_MAX_HIGHT，为了画图需要水位值赋为SEC_MAX_HIGHT
            re(n)=extraval;
        end   
    end
end