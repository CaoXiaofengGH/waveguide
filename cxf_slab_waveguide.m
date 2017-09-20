%************************导波光学-作业2-平板波导*********************************
%*********************************by曹晓峰 ********************************
%总体思路：
%1、导模特征方程左右式相减。在导模解附近的左右两侧，两个相减值异号。
%   两数相乘小于零即为异号。循环数值代入法找到较精确的异号点，其横坐标值即为导模的解
%2、为减少工作时长，采取方法是：
%   1、排除相减值较大的点（0.5）；
%   2、初步粗略循环判断两个相近（间隔0.01）的点是否异号。
%      找到初步异号点后，精确循环判断两个相近（0.0001）的点是否异号。
%说明：
%一
%功率章节中使用的场空间分布表达式与PPT前面章节的表达式不一致
%两者本质是一致的，处理方法上有差异，而主要参数是一致的，参数名替换后可以直接使用
%所以本程序算的是导波光学课本（曹庄琪）中的场空间分布
%即按照空间分段范围为(-inf,0),(0,d),(d,inf)
%
%二  
%本程序将各个参数的数据导入到了名为caoxiaofeng_waveguide.xlsx的EXCEL(2007)文件中
%可以在matlab的安装文件夹下bin里面找到
%一部分计算机可能不能运行数据导入部分程序，因为EXCEL的COM口被其他非系统程序占用
%但是程序其余部分依然可以正常运行。可以关闭EXCEL的com口的相关加载项来运行数据导出部分程序。
%如果感觉程序的数据导出部分影响很大，拖缓了程序的运行时间，或者影响图片的生成与显示
%可以将数据导出部分的语句注释掉，这不影响程序其他部分的运行
%因为这些数据也在命令行窗口输出。数据导出语句在程序的末尾处
%本程序在win7系统4G的RAM环境下15秒内可以运行完成

clear;
close all;
syms x X;
n1=2.0;             %芯层折射率
n2=1.5;             %衬底折射率
n3=1;               %覆盖层折射率
d=1.6e-6;           %芯层厚度
lambda=1.5e-6;      %光波波长
k=2*pi/lambda;      %波矢
c0=2.9979e8;        %真空光速
w0=k*c0;            %圆频率
u0=pi*4e-7;         %真空磁导率
eps0=8.8542e-12;    %真空介电常数

%导模特征方程左式
G(x)=tan(x);    

%TE模特征方程右式
V1=sqrt(n1^2-n2^2)*k*d;
V2=sqrt(n1^2-n3^2)*k*d;
F(x)=x*(sqrt(V1^2-x^2)+sqrt(V2^2-x^2))/(x^2-sqrt(V1^2-x^2)*sqrt(V2^2-x^2));

%TM模特征方程右式
U1=sqrt(n1^2-n3^2)*k*d;
U2=sqrt(n1^2-n2^2)*k*d;
H(x)=x*n1^2*(n2^2*sqrt(U1^2-x^2)+n3^2*sqrt(U2^2-x^2))/((n2*n3*x)^2-n1^4*sqrt(U1^2-x^2)*sqrt(U2^2-x^2));

%左右式相减
subE(x)=F(x)-G(x);
subM(x)=H(x)-G(x);

%***********特征方程左右式图形***********
xlenth=max(V1,U2)+0.5;
figure(1);
ezplot(G,[0,xlenth])
hold on;
ezplot(H,[0,xlenth])
hold on;
ezplot(F,[0,xlenth])
hold on;
title('Graphical');
xlabel('hd');


%******************求解TE导模个数********************
nmodeE=0;               %解的个数计数
solutionE=zeros(2,8);  %各个解，一般少于8个
for x=0.01:0.01:V1      %以间隔0.01初步循环判断
    if abs(subE(x))>0.5 %滤过左右式相差较大的点
        continue;
    end
    if eval(subE(x)*subE(x+0.01))<=0    %找到间隔为0.01的两个函数值异号的初步点
        tempe=x;          
        %从异号初步点开始以间隔0.0001精确循环判断，找到精确的异号点
        while eval(subE(tempe)*subE(tempe+0.0001))>0 
            tempe=tempe+0.0001;
        end
            solutionE(1,nmodeE+1)=tempe;     %将异号点的横坐标值即为导模的数值解
            solutionE(2,nmodeE+1)=tan(tempe);%异号点的纵坐标值
            nmodeE=nmodeE+1;               %计算解的个数
     end
end
solutionE

fprintf('共有%d个有效TE导模',nmodeE);              %表明有效解数量

%***********************求解TM导模个数******************
%*******************(注释与TE导模注释类似)***********
nmodeM=0;               %解的个数
solutionM=zeros(2,8);   %各个解
for x=0.01:0.01:U2
    if eval(abs(subM(x)))>0.5
        continue;
    end
    if eval(subM(x)*subM(x+0.01))<=0
        tempm=x;
        while eval(subM(tempm)*subM(tempm+0.0001))>0
            tempm=tempm+0.0001;
        end
         solutionM(1,nmodeM+1)=tempm;
         solutionM(2,nmodeM+1)=tan(tempm);
         nmodeM=nmodeM+1;
    end
end
solutionM

fprintf('共有%d个有效TM导模\n',nmodeM);

%*************图上标明解的位置****************
for j=1:nmodeE
    text(solutionE(1,j),solutionE(2,j),' E');
end
for j=1:nmodeM
    text(solutionM(1,j),solutionM(2,j),' M');
end

%***************TE各导模及其参数*****************
%不写成多维矩阵是为方便查看参数       
hE=zeros(1,nmodeE);         %h值
pE=zeros(1,nmodeE);         %p值
qE=zeros(1,nmodeE);         %q值
betaE=zeros(1,nmodeE);      %beta值
thetaE=zeros(1,nmodeE);     %theta值
heffE=zeros(1,nmodeE);      %heff值
AE=zeros(1,nmodeE);         %A值
phiE=zeros(1,nmodeE);       %场的空间分布表达式中cos(phi)中的phi值
sumPsubE=zeros(1,nmodeE);   %衬底层的总功率
sumPcoreE=zeros(1,nmodeE);  %芯层的总功率
sumPcoverE=zeros(1,nmodeE); %覆盖层的总功率

%求解TE模各参数及场的空间分布和能量的空间分布
%下面所计算的参数均是第三章PPT主体（前半部分内容）的表达式的各个参数
%下面所列的场的空间分布即是导波光学课本中的表达式
%即空间分布分段（-inf,0）,(0,d),(d,inf)
for j=1:nmodeE
    hE(j)=solutionE(1,j)/d;
    pE(j)=sqrt((n1^2-n2^2)*k^2-(solutionE(1,j)/d)^2);
    qE(j)=sqrt((n1^2-n3^2)*k^2-(solutionE(1,j)/d)^2);
    betaE(j)=sqrt(n1^2*k^2-(solutionE(1,j)/d)^2);
    thetaE(j)=asind(betaE(j)/(k*n1));
    heffE(j)=d+(1/pE(j))+(1/qE(j));
    AE(j)=2*hE(j)*sqrt(w0*u0/(heffE(j)*betaE(j)*(hE(j)^2+pE(j)^2)));
    phiE(j)=atan(pE(j)/hE(j));
    Esub=AE(j)*exp(pE(j)*X);                                            %TE模的衬底层电场空间分布
    Ecore=AE(j)*cos(hE(j)*X-phiE(j))/cos(phiE(j));                      %TE模的芯层电场空间分布
    Ecover=(AE(j)*cos(hE(j)*d-phiE(j))/cos(phiE(j)))*exp(-qE(j)*(X-d)); %TE模的覆盖层电场空间分布
    PsubE=betaE(j)*(Esub^2)/(2*w0*u0);                                  %TE模的衬底层电场空间分布
    PcoreE=betaE(j)*(Ecore^2)/(2*w0*u0);                                %TE模的芯层电场空间分布
    PcoverE=betaE(j)*(Ecover^2)/(2*w0*u0);                              %TE模的覆盖层电场空间分布
    
    %显示TE模的各层电场和能量的空间分布
    fprintf('第%d个TE导模\n',j);
    fprintf('Esub');
    vpa(Esub,5)
    fprintf('Ecore');
    vpa(Ecore,5)
    fprintf('Ecover');
    vpa(Ecover,5)
    
    fprintf('PsubE');
    vpa(PsubE,5)
    fprintf('PcoreE');
    vpa(PcoreE,5)
    fprintf('PcoverE');
    vpa(PcoverE,5)
end

%求解TE模各层的总功率
for j=1:nmodeE
    sumPsubE(j)=AE(j)^2*betaE(j)/(4*w0*u0*pE(j));
    sumPcoreE(j)=AE(j)^2*betaE(j)*(hE(j)^2+pE(j)^2)*(d+pE(j)/(hE(j)^2+pE(j)^2)+qE(j)/(hE(j)^2+qE(j)^2))/(4*w0*u0*hE(j)^2);
    sumPcoverE(j)=AE(j)^2*betaE(j)*(hE(j)^2+pE(j)^2)/(4*w0*u0*qE(j)*(hE(j)^2+qE(j)^2));
end

%显示TE模各参数
hE
pE
qE
betaE
thetaE
heffE
AE
sumPsubE
sumPcoreE
sumPcoverE

%显示TE模的电场和能量的空间分布
X1=-d:0.001*d:2*d;
for j=1:nmodeE
    Esub=AE(j)*exp(pE(j)*X1);
    Ecore=AE(j)*cos(hE(j)*X1-phiE(j))/cos(phiE(j));
    Ecover=(AE(j)*cos(hE(j)*d-phiE(j))/cos(phiE(j)))*exp(-qE(j)*(X1-d));
    PsubE=betaE(j)*power(Esub,2)/(2*w0*u0);
    PcoreE=betaE(j)*power(Ecore,2)/(2*w0*u0);
    PcoverE=betaE(j)*power(Ecover,2)/(2*w0*u0);
    PicE=Esub.*(X1<0)+Ecore.*((0<X1)&(X1<d))+Ecover.*(X1>d);        %电场的分段函数
    PicPE=PsubE.*(X1<0)+PcoreE.*((0<X1)&(X1<d))+PcoverE.*(X1>d);    %能量场的分段函数
    
    %分配各个画布，图形分割，标示
    figure(1+j);
    subplot(2,2,1);
    plot(X1,PicE)
    grid on;
    xlabel('x');
    ylabel('E');
    text(-d,AE(j)*exp(pE(j)*(-d)),'\leftarrow-d');
    text(0,AE(j)*exp(pE(j)*(0)),'\leftarrow0');
    text(d,(AE(j)*cos(hE(j)*d-phiE(j))/cos(phiE(j)))*exp(-qE(j)*(d-d)),'\leftarrowd');
    hold on;
    
    subplot(2,2,3);
    plot(X1,PicPE)
    grid on;
    xlabel('x');
    ylabel('P_E');
    text(-d,betaE(j)*power(AE(j)*exp(pE(j)*(-d)),2)/(2*w0*u0),'\leftarrow-d');
    text(0,betaE(j)*power(AE(j)*exp(pE(j)*0),2)/(2*w0*u0),'\leftarrow0');
    text(d,betaE(j)*power(AE(j)*cos(hE(j)*d-phiE(j))/cos(phiE(j)),2)/(2*w0*u0),'\leftarrowd');
    hold on;
end



%*******************TM各导模及其参数*****************
%*******************注释与TE模类似********************
hM=zeros(1,nmodeM);
pM=zeros(1,nmodeM);
qM=zeros(1,nmodeM);
betaM=zeros(1,nmodeM);
thetaM=zeros(1,nmodeM);
heffM=zeros(1,nmodeM);
AM=zeros(1,nmodeM);
phiM=zeros(1,nmodeM);
heffM=zeros(1,nmodeM);
sumPsubE=zeros(1,nmodeE);
sumPcoreE=zeros(1,nmodeE);
sumPcoverE=zeros(1,nmodeE);

%求解各参数
for j=1:nmodeM
    hM(j)=solutionM(1,j)/d;
    pM(j)=sqrt((n1^2-n2^2)*k^2-(solutionM(1,j)/d)^2);
    qM(j)=sqrt((n1^2-n3^2)*k^2-(solutionM(1,j)/d)^2);
    betaM(j)=sqrt(n1^2*k^2-(solutionM(1,j)/d)^2);
    thetaM(j)=asind(betaM(j)/(k*n1));
    heffM(j)=d+(n1^2*n3^2)*(hM(j)^2+qM(j)^2)/((n3^4*hM(j)^2+n1^4*qM(j)^2)*qM(j))+(n1^2*n2^2)*(hM(j)^2+pM(j)^2)/((n2^4*hM(j)^2+n1^4*pM(j)^2)*pM(j));
    AM(j)=2*hM(j)*sqrt((n1^2*n2^4)*w0*eps0/(betaM(j)*heffM(j)*(n2^4*hM(j)^2+n1^4*pM(j)^2)));
    phiM(j)=atan(n1^2*pM(j)/(n2^2*hM(j)));
    Hsub=AM(j)*exp(pM(j)*X);                                            %衬底层磁场函数
    Hcore=AM(j)*cos(hM(j)*X-phiM(j))/cos(phiM(j));                      %芯层磁场函数
    Hcover=(AM(j)*cos(hM(j)*d-phiM(j))/cos(phiM(j)))*exp(-qM(j)*(X-d)); %覆盖层磁场函数
    PsubM=betaM(j)*(Hsub^2)/(2*w0*eps0*n2^2);                           %衬底层能量场函数
    PcoreM=betaM(j)*(Hcore^2)/(2*w0*eps0*n1^2);                         %芯层能量场函数
    PcoverM=betaM(j)*(Hcover^2)/(2*w0*eps0*n3^2);                       %覆盖层能量场函数
    sumPsubM(j)=vpa(int(PsubM,X,-inf,0),5);                             %衬底层总能量
    sumPcoreM(j)=vpa(int(PcoreM,X,0,d),5);                              %芯层总能量
    sumPcoverM(j)=vpa(int(PcoverM,X,d,inf),5);                          %覆盖层总能量
    
    %显示磁场和能量场的空间分布函数
    fprintf('第%d个TM导模\n',j);
    fprintf('Hsub');
    vpa(Hsub,5)
    fprintf('PsubM');
    vpa(PsubM,5)
    fprintf('Hcore');
    vpa(Hcore,5)
    fprintf('PcoreM');
    vpa(PcoreM,5)
    fprintf('Hcover');
    vpa(Hcover,5)
    fprintf('PcoverM');
    vpa(PsubM,5)
end

%显示磁场与能量场的空间分布图形
X2=-d:0.001*d:2*d;
for j=1:3
    Hsub=AM(j)*exp(pM(j)*X2);
    Hcore=AM(j)*cos(hM(j)*X2-phiM(j))/cos(phiM(j));
    Hcover=(AM(j)*cos(hM(j)*d-phiM(j))/cos(phiM(j)))*exp(-qM(j)*(X2-d));
    PsubM=betaM(j)*power(Hsub,2)/(2*w0*eps0*n2^2);
    PcoreM=betaM(j)*power(Hcore,2)/(2*w0*eps0*n1^2);
    PcoverM=betaM(j)*power(Hcover,2)/(2*w0*eps0*n3^2);
    PicM=Hsub.*(X2<0)+Hcore.*(X2>0&X2<d)+Hcover.*(X2>d);        %磁场分布函数
    PicPM=PsubM.*(X2<0)+PcoreM.*(X2>0&X2<d)+PcoverM.*(X2>d);    %能量场分布函数
    
    %分配各个画布，图形分割，标示
    figure(1+j);
    subplot(2,2,2);
    plot(X2,PicM)
    grid on;
    xlabel('x');
    ylabel('H');
     
    text(-d,AM(j)*exp(pM(j)*(-d)),'\leftarrow-d');
    text(0,AM(j)*exp(pM(j)*0),'\leftarrow0');
    text(d,AM(j)*cos(hM(j)*d-phiM(j))/cos(phiM(j)),'\leftarrowd');
    hold on;
    
    subplot(2,2,4);
    plot(X2,PicPM)
    grid on;
    xlabel('x');
    ylabel('P_M');
    text(-d,betaM(j)*power(AM(j)*exp(pM(j)*(-d)),2)/(2*w0*eps0*n2^2),'\leftarrow-d');
    text(0,betaM(j)*power(AM(j)*exp(pM(j)*0),2)/(2*w0*eps0*n2^2),'\leftarrow0');
    text(d,betaM(j)*power(AM(j)*cos(hM(j)*d-phiM(j))/cos(phiM(j)),2)/(2*w0*eps0*n1^2),'\leftarrowd');
    hold on;
end

    
%*********显示TM模各参数**********

hM
pM
qM
betaM
thetaM
heffM
AM
sumPsubM
sumPcoreM
sumPcoverM


%**************************数据导出部分*************************
paranameE=['hE     ';'pE     ';'qE     ';'betaE  ';'thetaE ';'heffE  ';'AE     '];
paraE=zeros(7,nmodeE);
paranameM=['hM     ';'pM     ';'qM     ';'betaM  ';'thetaM ';'heffM  ';'AM     '];
paraM=zeros(7,nmodeM);

for j=1:nmodeE
    paraE(1,j)=hE(j);
    paraE(2,j)=pE(j);
    paraE(3,j)=qE(j);
    paraE(4,j)=betaE(j);
    paraE(5,j)=thetaE(j);
    paraE(6,j)=heffE(j);
    paraE(7,j)=AE(j);
end

for j=1:nmodeM
    paraM(1,j)=hM(j);
    paraM(2,j)=pM(j);
    paraM(3,j)=qM(j);
    paraM(4,j)=betaM(j);
    paraM(5,j)=thetaM(j);
    paraM(6,j)=heffM(j);
    paraM(7,j)=AM(j);
end

xlswrite('caoxiaofeng_waveguide.xlsx',paranameE,'sheet1','A1')
xlswrite('caoxiaofeng_waveguide.xlsx',paraE,'sheet1','I1')
xlswrite('caoxiaofeng_waveguide.xlsx',paranameM,'sheet1','A9')
xlswrite('caoxiaofeng_waveguide.xlsx',paraM,'sheet1','I9')
%******************************数据导出完毕**********************
