%*******************HOMEWORK4*******************
%****************计算光纤EH模的截止频率**********
clear;
close all;

a=10e-6;    %纤芯半径
n1=1.49;    %纤芯折射率
n2=1.489;   %薄层折射率

%**********************数值求解*****************
besselj1=@(x)besselj(1,x);
for n=1:5
    Vc1(n)=fzero(besselj1,[(n-1) n]*pi);    %归一化截止频率
    lambda_cut(n)=2*pi*a*sqrt(n1^2-n2^2)/Vc1(n);    %截止波长
end

%**********bessel函数图和零点********************
x=0:pi/100:10*pi;
y=besselj(1,x);
plot(Vc1,zeros(1,5),'o',x,y)
line([0 6*pi],[0 0],'color','black');
axis([0 6*pi -0.5 1.0]);
xlabel('U');
ylabel('J_1(U)');
title('一阶贝塞尔函数');
text(Vc1(3),-0.1,'j_1_,_3');

Vc1
lambda_cut
