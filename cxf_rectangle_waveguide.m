%*************************导波光学-矩形波导_Homework3**********************
%*********************************by曹晓峰 ********************************
clear;
close all;

syms Kx Ky X Y;

%初始参数
n1=1.4549;  %矩形中心折射率
n2=1.444;   %上侧折射率
n3=1.444;   %下侧
n4=1.444;   %左侧
n5=1.444;   %右侧

a=4e-6;         %矩形长度
b=4e-6;         %矩形宽度
lambda=1.5e-6;  %波长
k0=2*pi/lambda; %波数
c0=2.99792e8;   %真空光速
w=k0*c0;        %圆频率
u0=pi*4e-7;     %这空介电常数

%衰减因子
alpha4x(Kx)=sqrt((k0^2)*(n1^2-n4^2)-Kx^2);
alpha5x(Kx)=sqrt((k0^2)*(n1^2-n5^2)-Kx^2);
alpha2y(Ky)=sqrt((k0^2)*(n1^2-n2^2)-Ky^2);
alpha3y(Ky)=sqrt((k0^2)*(n1^2-n3^2)-Ky^2);

%K值的有效取值范围
limitKx=min(k0*sqrt(n1^2-n4^2),k0*sqrt(n1^2-n5^2));
limitKy=min(k0*sqrt(n1^2-n2^2),k0*sqrt(n1^2-n3^2));

%色散方程
G(Kx)=tan(2*Kx*a);
ExF(Kx)=n1^2*Kx*(alpha4x(Kx)*n5^2+alpha5x(Kx)*n4^2)/((n4*n5*Kx)^2-(n1^4*alpha4x(Kx)*alpha5x(Kx)));
R(Ky)=tan(2*Ky*b);
ExS(Ky)=Ky*(alpha2y(Ky)+alpha3y(Ky))/(Ky^2-alpha2y(Ky)*alpha3y(Ky));

%方程左右两式相减
subxEx(Kx)=ExF(Kx)-G(Kx);
subyEx(Ky)=ExS(Ky)-R(Ky);

%画色散方程图形
figure(1);
ezplot(G(Kx),[0,1.05*limitKx])
hold on;
ezplot(ExF(Kx),[0,1.05*limitKx])
hold on;
title('Graphic');

figure(2);
ezplot(R(Ky),[0,1.05*limitKy])
hold on;
ezplot(ExS(Ky),[0,1.05*limitKy])
hold on;
title('Graphic');

%******************求解Ex模式的本征方程中的Kx********************
nKx=0;                                         %解的个数计数
solutionKx=zeros(2,8);                         %各个解,一般少于8个
for Kx=0.01*limitKx:0.01*limitKx:limitKx       %以间隔0.01初步循环判断
    if abs(subxEx(Kx))>2                       %滤过左右式相差较大的点
        continue;
    end
    if eval(subxEx(Kx)*subxEx(Kx+0.01*limitKx))<=0    %找到间隔为0.01的两个函数值异号的初步点
        tempe=Kx;
           %从异号初步点开始以间隔0.0001精确循环判断
        while eval(subxEx(tempe)*subxEx(tempe+0.0001*limitKx))>0    %循环判断找到精确的异号点
            tempe=tempe+0.0001*limitKx;
        end
            solutionKx(1,nKx+1)=tempe;          %将异号点的横坐标值即为导模的数值解
            solutionKx(2,nKx+1)=tan(2*a*tempe);   %异号点的纵坐标值
            nKx=nKx+1;                          %计算解的个数
    end
end
solutionKx                                      %显示各个解
fprintf('Ex模式共有%d个有效Kx值\n',nKx);         %表明有效解数量

%******************求解Ex模式的本征方程中的Ky********************
nKy=0;                                          %解的个数计数
solutionKy=zeros(2,8);                          %各个解,一般少于8个
for Ky=0.01*limitKy:0.01*limitKy:limitKy        %以间隔0.01初步循环判断
    if abs(subyEx(Ky))>2                        %滤过左右式相差较大的点
        continue;
    end
    if eval(subyEx(Ky)*subyEx(Ky+0.01*limitKy))<=0    %找到间隔为0.01的两个函数值异号的初步点
        tempe=Ky;
           %从异号初步点开始以间隔0.0001精确循环判断
        while eval(subyEx(tempe)*subyEx(tempe+0.0001*limitKy))>0    %循环判断找到精确的异号点
            tempe=tempe+0.0001*limitKy;
        end
            solutionKy(1,nKy+1)=tempe;        %将异号点的横坐标值即为导模的数值解
            solutionKy(2,nKy+1)=tan(2*b*tempe);   %异号点的纵坐标值
            nKy=nKy+1;                    %计算解的个数
    end
end
solutionKy                                         %显示各个解
fprintf('Ex模式共有%d个有效Ky值\n',nKy);              %表明有效解数量

%*******************计算衰减系数/相位/幅度************************
v_alpha4x=zeros(1,nKx);             %衰减系数
v_alpha5x=zeros(1,nKx);
v_alpha2y=zeros(1,nKy);
v_alpha3y=zeros(1,nKy);
epsi=zeros(1,nKx);                  %初始相位
eta=zeros(1,nKy);
A1=1;                               %磁场振幅因子
A4=zeros(1,nKx);
A5=zeros(1,nKx);
A2=zeros(1,nKy);
A3=zeros(1,nKy);

for j=1:nKx                         %alpha4x,alpha5x的值
    v_alpha4x(j)=alpha4x(solutionKx(1,j));
    v_alpha5x(j)=alpha5x(solutionKx(1,j));
    epsi(j)=-solutionKx(1,j)*a+atan((n1^2*v_alpha4x(j))/(n4^2*solutionKx(1,j)));
    A4(j)=A1*cos(a*solutionKx(1,j)+epsi(j))/(exp(-a*v_alpha4x(j)));
    A5(j)=A1*cos(a*solutionKx(1,j)-epsi(j))/(exp(-a*v_alpha5x(j)));
end

for t=1:nKy                         %alpha2y,alpha3y的值
    v_alpha2y(t)=alpha2y(solutionKy(1,t));
    v_alpha3y(t)=alpha3y(solutionKy(1,t));
    eta(t)=-solutionKy(1,t)*b+atan(v_alpha3y(t)/solutionKy(1,t));
    A2(t)=A1*cos(b*solutionKy(1,t)-eta(t))/exp(-b*v_alpha2y(t));
    A3(t)=A1*cos(b*solutionKy(1,t)+eta(t))/exp(-b*v_alpha3y(t));
end

%**************************各个模式下磁场的分布*********************
beta=zeros(nKx,nKy);
for j=1:nKx
    for t=1:nKy
        %传播常数
        beta(j,t)=sqrt((n1*k0)^2-solutionKx(j)^2-solutionKy(t)^2);
        
        %y方向磁场
        H1(X,Y)=A1*cos(solutionKx(1,j)*X+epsi(j))*cos(solutionKy(1,t)*Y+eta(t));
        H2(X,Y)=A2(t)*cos(solutionKx(1,j)*X+epsi(j))*exp(v_alpha2y(t)*Y);
        H3(X,Y)=A3(t)*cos(solutionKx(1,j)*X+epsi(j))*exp(-v_alpha3y(t)*Y);
        H4(X,Y)=A4(j)*cos(solutionKy(1,t)*Y+eta(t))*exp(-v_alpha4x(j)*X);
        H5(X,Y)=A5(j)*cos(solutionKy(1,t)*Y+eta(t))*exp(v_alpha5x(j)*X);
        
        %x方向电场
        E1(X,Y)=H1(X,Y)*(w*u0-solutionKx(j)^2/(w*u0*(n1^2)))/beta(j,t);
        E2(X,Y)=H2(X,Y)*(w*u0-solutionKx(j)^2/(w*u0*(n2^2)))/beta(j,t);
        E3(X,Y)=H3(X,Y)*(w*u0-solutionKx(j)^2/(w*u0*(n3^2)))/beta(j,t);
        E4(X,Y)=H4(X,Y)*(w*u0+alpha4x(j)^2/(w*u0*(n4^2)))/beta(j,t);
        E5(X,Y)=H5(X,Y)*(w*u0+alpha5x(j)^2/(w*u0*(n5^2)))/beta(j,t);
        
        %显示磁场分布函数
        fprintf('Ex的%d%d模式\n',j,t);
        fprintf('H1\n');
        vpa(H1,5)
        fprintf('H2\n');
        vpa(H2,5)
        fprintf('H3\n');
        vpa(H3,5)
        fprintf('H4\n');
        vpa(H4,5)
        fprintf('H5\n');
        vpa(H5,5)
        
        %显示电场分布函数
        fprintf('E1\n');
        vpa(E1,5)
        fprintf('E2\n');
        vpa(E2,5)
        fprintf('E3\n');
        vpa(E3,5)
        fprintf('E4\n');
        vpa(E4,5)
        fprintf('E5\n');
        vpa(E5,5)
        
        %画出各个模式下磁场的分布
        figure(j+2*t);                      %分画布
        x=-2*a:0.02*a:2*a;
        y=-2*b:0.02*b:2*b;
        [xx,yy]=meshgrid(x,y);              %建立网格
        Pic_H1=(A1*cos(solutionKx(1,j)*xx+epsi(j)).*cos(solutionKy(1,t)*yy+eta(t)));
        Pic_H2=(A2(t)*cos(solutionKx(1,j)*xx+epsi(j)).*exp(v_alpha2y(t)*yy));
        Pic_H3=(A3(t)*cos(solutionKx(1,j)*xx+epsi(j)).*exp(-v_alpha3y(t)*yy));
        Pic_H4=(A4(j)*cos(solutionKy(1,t)*yy+eta(t)).*exp(-v_alpha4x(j)*xx));
        Pic_H5=(A5(j)*cos(solutionKy(1,t)*yy+eta(t)).*exp(v_alpha5x(j)*xx));
        Pic_H=Pic_H1.*((xx>=-a)&(xx<=a)&(yy>=-b)&(yy<=b))+Pic_H2.*((xx>-a)&(xx<a)&(yy<-b))+Pic_H3.*((xx>-a)&(xx<a)&(yy>b))+Pic_H4.*((yy>-b)&(yy<b)&(xx>a))+Pic_H5.*((yy>-b)&(yy<b)&(xx<-a));
        meshc(xx,yy,Pic_H)                  %画三维网线图 
        hidden off;
        hold on;
        title('H');
        xlabel('x');
        ylabel('y');
        
    figure(nKx*nKy+j+2*t); 
        u=-2*a:0.02*a:2*a;
        v=-2*b:0.02*b:2*b;
        [uu,vv]=meshgrid(u,v);    
        Pic_E1=((w*u0-solutionKx(j)^2/(w*u0*(n1^2)))/beta(j,t))*(A1*cos(solutionKx(1,j)*uu+epsi(j)).*cos(solutionKy(1,t)*vv+eta(t)));
        Pic_E2=((w*u0-solutionKx(j)^2/(w*u0*(n2^2)))/beta(j,t))*(A2(t)*cos(solutionKx(1,j)*uu+epsi(j)).*exp(v_alpha2y(t)*vv));
        Pic_E3=((w*u0-solutionKx(j)^2/(w*u0*(n3^2)))/beta(j,t))*(A3(t)*cos(solutionKx(1,j)*uu+epsi(j)).*exp(-v_alpha3y(t)*vv));
        Pic_E4=((w*u0+v_alpha4x(j)^2/(w*u0*(n4^2)))/beta(j,t))*(A4(j)*cos(solutionKy(1,t)*vv+eta(t)).*exp(-v_alpha4x(j)*uu));
        Pic_E5=((w*u0+v_alpha5x(j)^2/(w*u0*(n5^2)))/beta(j,t))*(A5(j)*cos(solutionKy(1,t)*vv+eta(t)).*exp(v_alpha5x(j)*uu));
        Pic_E=Pic_E1.*((uu>=-a)&(uu<=a)&(vv>=-b)&(vv<=b))+Pic_E2.*((uu>-a)&(uu<a)&(vv<-b))+Pic_E3.*((uu>-a)&(uu<a)&(vv>b))+Pic_E4.*((vv>-b)&(vv<b)&(uu>a))+Pic_E5.*((vv>-b)&(vv<b)&(uu<-a));
        meshc(uu,vv,Pic_E)
        hidden off;
        hold on;
        title('E');
        xlabel('x');
        ylabel('y');
        
    end
end

beta
epsi
eta
