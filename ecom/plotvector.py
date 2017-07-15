function [ArrowID,AllArrowMaxX]=plotvector(X,Y,Vel,Dir,userscale,arrowlength,ColorStr,FilltagStr,UVTagStr,ArrowLenModeStr,MaxVectorLen,AxisID)
%--------------------------------------------------------------------------
%[ArrowID]=plotvector(X,Y,Vel,Dir,userscale,arrowlength,ColorStr,FilltagStr,UVTagStr,ArrowLenModeStr,MaxVectorLen,AxisID)
%--------------------------------------------------------------------------
%    绘制流矢图，功能同quiver命令，但强于quiver命令
%              李为华, 2011.10.16
%             liweihua0903@163.com
%--------------------------------------------------------------------------
%优点：
%1.自动根据当前绘图画布x、y实际比例调整dx、dy，当x、y轴绘图比尺不同时，可避免流矢长度及方向失真问题，
%   尤为适用于绘制随时间变化的定点流速测量结果
%2.仅输出当前绘图画布中x、y限定范围内的流矢，以减小绘图量及成图存储空间
%3.可通过箭头长度设定方式，设定是否绘制矢量箭头及箭头长度大小（流速越大，箭头越长，最长为平均流速对应箭头）
%4.可采用填充三角形方式绘制矢量箭头
%5.可方便绘制流场叠加效果图
%-----------------------------------------------
%变量说明：
% X、Y、Vel、Dir分别为x、y坐标和流速、流向（正北为零，顺时针方向）
% userscale,流矢长度比尺，相对参数，数值越大，流矢长度越大
% arrowlength,箭头长度设定参数，相对参数，数值越大，箭头相对越长
% ColorStr,流矢颜色，字符串变量，只能为单色设置，如 'k'
% FilltagStr,是否采用填充三角形式绘制箭头标志，取值为'fill'时采用三角填充，取值为'line'时箭头为线形
%            当绘制大面积流场使用填充模式绘图时，计算机可能会有一定响应延时，请耐心等待
% UVTagStr,取值为"uv"时，第三、第四变量为u、v分量，取值为“vd”为流速流向
% ArrowLenModeStr，取值为"fix"时，所有箭头长度等长，取值为"var"时箭头长度随矢量大小而变化
% MaxVectorLen，最大矢量长度，如设定最大流速为5(m/s)，则相应sqrt(x宽度^2+y宽度^2)/80*userscale代表5m/s;
%             如由程序自动决定矢量长度，则该变量可不输入或设定为一不大于0的数值；
%         注意：
%             (1)当需要进行流场叠加显示时，该参数必须设定为可靠地正数值才能保证前后比尺相同；
%             (2)当该变量设定数值小于max(Vel)时，长于该设定数值的矢量将限定为该数值长度。
% AxisID,绘图坐标轴ID
%-----------------------------------------------
% 调用方式：
% 直接使用流速流向时,
% plotvector(X,Y,Vel,Dir)
% plotvector(X,Y,Vel,Dir,userscale)
% plotvector(X,Y,Vel,Dir,userscale,arrowlength)
% plotvector(X,Y,Vel,Dir,userscale,arrowlength,ColorStr)
% plotvector(X,Y,Vel,Dir,userscale,arrowlength,ColorStr,FilltagStr)
% plotvector(X,Y,Vel,Dir,userscale,arrowlength,ColorStr,FilltagStr,UVTagStr)
% plotvector(X,Y,Vel,Dir,userscale,arrowlength,ColorStr,FilltagStr,UVTagStr,ArrowLenModeStr)
% plotvector(X,Y,Vel,Dir,userscale,arrowlength,ColorStr,FilltagStr,UVTagStr,ArrowLenModeStr,MaxVectorLen)
% 直接使用U、V分量时，
% plotvector(X,Y,dx,dy,userscale,arrowlength,ColorStr,FilltagStr,'UV')
% plotvector(X,Y,dx,dy,userscale,arrowlength,ColorStr,FilltagStr,'UV',ArrowLenModeStr)
% plotvector(X,Y,dx,dy,userscale,arrowlength,ColorStr,FilltagStr,'UV',ArrowLenModeStr,MaxVectorLen)
% plotvector(X,Y,dx,dy,userscale,arrowlength,ColorStr,FilltagStr,'UV',ArrowLenModeStr,MaxVectorLen,AxisID)
% 如果userscale=1，则坐标轴宽度和高度最小值的10%长度作为最大矢量长度，缺省默认值
% 如果arrowlength=1,则坐标轴宽度和高度最小值的2%长度作为最大箭头长度，缺省默认值
% 程序缺省默认值：
%    userscale=1;
%    arrowlength=1;
%    ColorStr='b';
%    FilltagStr='line';
%    UVTagStr='vd';
%    ArrowLenModeStr='fix';
%    MaxVectorLen=-1;
%    AxisID=gca;
%--------------------------------------------------------------------------
%Example:
% [x,y] = meshgrid(-2:.2:2,-1:.15:1);
% z = x .* exp(-x.^2 - y.^2); [px,py] = gradient(z,.2,.15);
% figure(1);contour(x,y,z); hold on;
% a=get(gcf,'position');a(1)=a(1)-40;a(2)=a(2)-a(4);a(4)=a(4)*2;
% set(gcf,'position',a);
% quiver(x,y,px,py); set(gca,'xlim',[-2 2],'ylim',[-1 1]);
% title('quiver(蓝色)和plotvector(红色)效果比较');
% plotvector(x,y,px,py,1,1,'r','fill','uv','fix',1); hold off;
% figure(3);contour(x,y,z); hold on;
% a=get(gcf,'position');a(2)=a(2)-a(4);a(3)=a(3)*2;
% set(gcf,'position',a);
% quiver(x,y,px,py);axis equal;
% plotvector(x,y,px,py,0.8,1,'r','line','uv','fix',1); hold off;
% title('两命令正常显示效果叠加，quiver蓝色，plotvector红色');
%--------------------------------------------------------------------------
Angle = pi / 18;  %箭头与流矢夹角预设为10°
%--------------------------------------------------------------------------
%输入变量处理
%--------------------------------------------------------------------------
if nargin<4
    error('变量输入信息不足');
elseif nargin==4
    arrowlength=1;
    userscale=1;
    ColorStr='b';
    FilltagStr='line';
    UVTagStr='vd';
    ArrowLenModeStr='fix';
    MaxVectorLen=-1;
    AxisID=gca;
elseif nargin==5
    arrowlength=1;
    ColorStr='b';
    FilltagStr='line';
    UVTagStr='vd';
    ArrowLenModeStr='fix';
    MaxVectorLen=-1;
    AxisID=gca;
elseif nargin==6
    ColorStr='b';
    FilltagStr='line';
    UVTagStr='vd';
    ArrowLenModeStr='fix';
    MaxVectorLen=-1;
    AxisID=gca;
elseif nargin==7
    FilltagStr='line';
    UVTagStr='vd';
    ArrowLenModeStr='fix';
    MaxVectorLen=-1;
    AxisID=gca;
elseif nargin==8
    UVTagStr='vd';
    ArrowLenModeStr='fix';
    MaxVectorLen=-1;
    AxisID=gca;
elseif nargin==9
    ArrowLenModeStr='fix';
    MaxVectorLen=-1;
    AxisID=gca;
elseif nargin==10
    MaxVectorLen=-1;
    AxisID=gca;
elseif nargin==11
    AxisID=gca;
end
axes(AxisID);
Filltag=strcmpi(FilltagStr,'fill');
UVTag=strcmpi(UVTagStr,'uv');
ArrowLenMode=strcmpi(ArrowLenModeStr,'fix');

if UVTag==1
    dx=Vel;
    dy=Dir;
    [Vel,Dir]=EN2VD(dx,dy);
else
    dx=Vel.*sin(Dir.*pi./180);
    dy=Vel.*cos(Dir.*pi./180);
end
%--------------------------------------------------------------------------
%画布处理，保持绘图前的宽度、高度、宽高比、xy轴显示范围和hold状态
%--------------------------------------------------------------------------
tempvar1=get(gcf,'position');
tempvar2=get(AxisID,'position');
XAxisLimit=get(AxisID,'xlim');
YAxisLimit=get(AxisID,'ylim');
plotstatus=ishold(AxisID); %当前hold状态
if (XAxisLimit(1)==0)&(YAxisLimit(1)==0)&(XAxisLimit(2)==1)&(YAxisLimit(2)==1)
    IsNewPlot=1;
else
    IsNewPlot=0;
end
if IsNewPlot==1
    if plotstatus==0
        TempID=quiver(AxisID,X,Y,dx,dy);
        %if Filltag==1
            %Filltag=0;
            %warning('程序自动新创建绘图窗体时不可使用fill模式绘制箭头。');
        %end
    end
end
tempvar1=get(gcf,'position');
tempvar2=get(AxisID,'position');
XAxisLimit=get(AxisID,'xlim');
YAxisLimit=get(AxisID,'ylim');
xwidth=tempvar1(3)*tempvar2(3); %画布宽度
ywidth=tempvar1(4)*tempvar2(4); %画布高度
tempvar3=xwidth/(XAxisLimit(2)-XAxisLimit(1));
tempvar4=ywidth/(YAxisLimit(2)-YAxisLimit(1));
xscale=tempvar3/(tempvar3+tempvar4);
yscale=tempvar4/(tempvar3+tempvar4);
if IsNewPlot==1
    if plotstatus==0
        delete(TempID);
        plot(AxisID,nan,nan);
        set(AxisID,'xlim',XAxisLimit,'ylim',YAxisLimit);
    end
end
%--------------------------------------------------------------------------
%绘图数据处理
%--------------------------------------------------------------------------
if (IsNewPlot==0) %剔除不在当前x、y轴显示范围内的数据点
    Tag_Valid=find((X>=XAxisLimit(1))&(X<=XAxisLimit(2))&(Y>=YAxisLimit(1))&(Y<=YAxisLimit(2)));
    X=X(Tag_Valid);
    Y=Y(Tag_Valid);
    Vel=Vel(Tag_Valid);
    Dir=Dir(Tag_Valid);
end
%u、v数据根据当前画布x、y比例尺相应变换，以纠正xy比例尺不同导致的矢量变形
%并调整矢量长度至设定长度
BaseVectorSize=min([(XAxisLimit(2)-XAxisLimit(1)),(YAxisLimit(2)-YAxisLimit(1))])./10; %坐标轴宽度和高度最小值的10%长度作为最大矢量长度
n=length(X);

if MaxVectorLen>0 %设定超过限制的矢量长度为限制长度
    TempVar1=find(Vel>MaxVectorLen);
    if isempty(TempVar1)~=1
        Vel(TempVar1)=MaxVectorLen;
    end
end
Vel=Vel.*BaseVectorSize./max([Vel(:);MaxVectorLen]) .* userscale;  %调整矢量单位至绘图单位
dx=Vel.*sin(Dir.*pi./180);
dy=Vel.*cos(Dir.*pi./180);
%箭头长度处理
BaseArrowSize=BaseVectorSize./10;%最大矢量长度的10%作为最大箭头长度
if ArrowLenMode==1 %箭头大小固定
    L_Arrow=BaseArrowSize.*arrowlength;
end
%箭头填充方式预处理
if Filltag==0
    ArrowX=zeros(7*n,1);
    ArrowY=zeros(7*n,1);
else
    ArrowX=zeros(3*n,1);
    ArrowY=zeros(3*n,1);
    FillX=zeros(3,1);
    FillY=zeros(3,1);
end
%生成流矢绘图数据

AllArrowMaxX=dx.*0;

for i=1:n
    cosA = dx(i) / sqrt(dx(i) ^ 2 + dy(i) ^ 2);
    sinA = dy(i) / sqrt(dx(i) ^ 2 + dy(i) ^ 2);
    if ArrowLenMode==0 %箭头大小随流矢长度变化而变化
        L_Arrow = Vel(i) * arrowlength;
    end
    if (L_Arrow>Vel(i))&(ArrowLenMode==1)
        TempLen_Arrow=Vel(i).*0.7;
    else
        TempLen_Arrow=L_Arrow;
    end
    arowx_up =  TempLen_Arrow * (cosA * cos(-Angle) - sinA * sin(-Angle));
    arowy_up =  TempLen_Arrow * (sinA * cos(-Angle) + cosA * sin(-Angle));
    arowx_down = TempLen_Arrow * (cosA * cos(Angle) - sinA * sin(Angle));
    arowy_down = TempLen_Arrow * (sinA * cos(Angle) + cosA * sin(Angle));
    if Filltag==0
        ArrowX((i-1)*7+1,1)=X(i);
        ArrowY((i-1)*7+1,1)=Y(i);
        ArrowX((i-1)*7+2,1)=X(i)+dx(i)./ xscale;
        ArrowY((i-1)*7+2,1)=Y(i)+dy(i)./ yscale;
        ArrowX((i-1)*7+3,1)=nan;
        ArrowY((i-1)*7+3,1)=nan;
        ArrowX((i-1)*7+4,1)=X(i)+(dx(i) - arowx_up)./ xscale;
        ArrowY((i-1)*7+4,1)=Y(i)+(dy(i) - arowy_up)./ yscale;
        ArrowX((i-1)*7+5,1)=X(i)+dx(i)./ xscale;
        ArrowY((i-1)*7+5,1)=Y(i)+dy(i)./ yscale;
        ArrowX((i-1)*7+6,1)=X(i)+(dx(i) - arowx_down)./ xscale;
        ArrowY((i-1)*7+6,1)=Y(i)+(dy(i) - arowy_down)./ yscale;
        ArrowX((i-1)*7+7,1)=nan;
        ArrowY((i-1)*7+7,1)=nan;
    else
        ArrowX((i-1)*3+1,1)=X(i);
        ArrowY((i-1)*3+1,1)=Y(i);
        ArrowX((i-1)*3+2,1)=X(i)+dx(i)./ xscale;
        ArrowY((i-1)*3+2,1)=Y(i)+dy(i)./ yscale;
        ArrowX((i-1)*3+3,1)=nan;
        ArrowY((i-1)*3+3,1)=nan;
        
        FillX(1,1)=X(i)+(dx(i) - arowx_up)./ xscale;
        FillY(1,1)=Y(i)+(dy(i) - arowy_up)./ yscale;
        FillX(2,1)=X(i)+dx(i)./ xscale;
        FillY(2,1)=Y(i)+dy(i)./ yscale;
        FillX(3,1)=X(i)+(dx(i) - arowx_down)./ xscale;
        FillY(3,1)=Y(i)+(dy(i) - arowy_down)./ yscale;
        if (plotstatus==0)&&(i==2)
            hold(AxisID,'on');
            set(AxisID,'xlim',XAxisLimit,'ylim',YAxisLimit);
        end
        FillID=fill(FillX,FillY,ColorStr);
        set(FillID,'edgecolor',ColorStr);
    end
    AllArrowMaxX(i)=X(i)+dx(i)./ xscale;
end
%--------------------------------------------------------------------------
%绘图
%--------------------------------------------------------------------------
ArrowID=plot(AxisID,ArrowX,ArrowY,ColorStr);
set(AxisID,'xlim',XAxisLimit,'ylim',YAxisLimit);

if plotstatus==0
    hold(AxisID,'off'); 
else
    hold(AxisID,'on'); 
end
% end of plotvector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [V,D]=EN2VD(V_E,V_N)
%由流速的东/北分量计算流速流向
V=sqrt(V_E.^2+V_N.^2);
c2=find((V_E>=0)&(V_N>=0));  %第一象限内流向
D(c2)=asin(abs(V_E(c2))./V(c2)).*180./pi;
clear c2;
c2=find((V_E>=0)&(V_N<0));  %第四象限内流向
D(c2)=180-asin(abs(V_E(c2))./V(c2)).*180./pi;
clear c2;
c2=find((V_E<0)&(V_N>=0));  %第二象限内流向
D(c2)=360-asin(abs(V_E(c2))./V(c2)).*180./pi;
clear c2;
c2=find((V_E<0)&(V_N<0));  %第三象限内流向
D(c2)=180+asin(abs(V_E(c2))./V(c2)).*180./pi;
clear c2;
D=D';