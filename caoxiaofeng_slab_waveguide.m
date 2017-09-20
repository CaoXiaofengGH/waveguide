%************************������ѧ-��ҵ2-ƽ�岨��*********************************
%*********************************by������ ********************************
%����˼·��
%1����ģ������������ʽ������ڵ�ģ�⸽�����������࣬�������ֵ��š�
%   �������С���㼴Ϊ��š�ѭ����ֵ���뷨�ҵ��Ͼ�ȷ����ŵ㣬�������ֵ��Ϊ��ģ�Ľ�
%2��Ϊ���ٹ���ʱ������ȡ�����ǣ�
%   1���ų����ֵ�ϴ�ĵ㣨0.5����
%   2����������ѭ���ж�������������0.01���ĵ��Ƿ���š�
%      �ҵ�������ŵ�󣬾�ȷѭ���ж����������0.0001���ĵ��Ƿ���š�
%˵����
%һ
%�����½���ʹ�õĳ��ռ�ֲ����ʽ��PPTǰ���½ڵı��ʽ��һ��
%���߱�����һ�µģ����������в��죬����Ҫ������һ�µģ��������滻�����ֱ��ʹ��
%���Ա���������ǵ�����ѧ�α�����ׯ�����еĳ��ռ�ֲ�
%�����տռ�ֶη�ΧΪ(-inf,0),(0,d),(d,inf)
%
%��  
%�����򽫸������������ݵ��뵽����Ϊcaoxiaofeng_waveguide.xlsx��EXCEL(2007)�ļ���
%������matlab�İ�װ�ļ�����bin�����ҵ�
%һ���ּ�������ܲ����������ݵ��벿�ֳ�����ΪEXCEL��COM�ڱ�������ϵͳ����ռ��
%���ǳ������ಿ����Ȼ�����������С����Թر�EXCEL��com�ڵ���ؼ��������������ݵ������ֳ���
%����о���������ݵ�������Ӱ��ܴ��ϻ��˳��������ʱ�䣬����Ӱ��ͼƬ����������ʾ
%���Խ����ݵ������ֵ����ע�͵����ⲻӰ������������ֵ�����
%��Ϊ��Щ����Ҳ�������д�����������ݵ�������ڳ����ĩβ��
%��������win7ϵͳ4G��RAM������15���ڿ����������

clear;
close all;
syms x X;
n1=2.0;             %о��������
n2=1.5;             %�ĵ�������
n3=1;               %���ǲ�������
d=1.6e-6;           %о����
lambda=1.5e-6;      %�Ⲩ����
k=2*pi/lambda;      %��ʸ
c0=2.9979e8;        %��չ���
w0=k*c0;            %ԲƵ��
u0=pi*4e-7;         %��մŵ���
eps0=8.8542e-12;    %��ս�糣��

%��ģ����������ʽ
G(x)=tan(x);    

%TEģ����������ʽ
V1=sqrt(n1^2-n2^2)*k*d;
V2=sqrt(n1^2-n3^2)*k*d;
F(x)=x*(sqrt(V1^2-x^2)+sqrt(V2^2-x^2))/(x^2-sqrt(V1^2-x^2)*sqrt(V2^2-x^2));

%TMģ����������ʽ
U1=sqrt(n1^2-n3^2)*k*d;
U2=sqrt(n1^2-n2^2)*k*d;
H(x)=x*n1^2*(n2^2*sqrt(U1^2-x^2)+n3^2*sqrt(U2^2-x^2))/((n2*n3*x)^2-n1^4*sqrt(U1^2-x^2)*sqrt(U2^2-x^2));

%����ʽ���
subE(x)=F(x)-G(x);
subM(x)=H(x)-G(x);

%***********������������ʽͼ��***********
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


%******************���TE��ģ����********************
nmodeE=0;               %��ĸ�������
solutionE=zeros(2,8);  %�����⣬һ������8��
for x=0.01:0.01:V1      %�Լ��0.01����ѭ���ж�
    if abs(subE(x))>0.5 %�˹�����ʽ���ϴ�ĵ�
        continue;
    end
    if eval(subE(x)*subE(x+0.01))<=0    %�ҵ����Ϊ0.01����������ֵ��ŵĳ�����
        tempe=x;          
        %����ų����㿪ʼ�Լ��0.0001��ȷѭ���жϣ��ҵ���ȷ����ŵ�
        while eval(subE(tempe)*subE(tempe+0.0001))>0 
            tempe=tempe+0.0001;
        end
            solutionE(1,nmodeE+1)=tempe;     %����ŵ�ĺ�����ֵ��Ϊ��ģ����ֵ��
            solutionE(2,nmodeE+1)=tan(tempe);%��ŵ��������ֵ
            nmodeE=nmodeE+1;               %�����ĸ���
     end
end
solutionE

fprintf('����%d����ЧTE��ģ',nmodeE);              %������Ч������

%***********************���TM��ģ����******************
%*******************(ע����TE��ģע������)***********
nmodeM=0;               %��ĸ���
solutionM=zeros(2,8);   %������
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

fprintf('����%d����ЧTM��ģ\n',nmodeM);

%*************ͼ�ϱ������λ��****************
for j=1:nmodeE
    text(solutionE(1,j),solutionE(2,j),' E');
end
for j=1:nmodeM
    text(solutionM(1,j),solutionM(2,j),' M');
end

%***************TE����ģ�������*****************
%��д�ɶ�ά������Ϊ����鿴����       
hE=zeros(1,nmodeE);         %hֵ
pE=zeros(1,nmodeE);         %pֵ
qE=zeros(1,nmodeE);         %qֵ
betaE=zeros(1,nmodeE);      %betaֵ
thetaE=zeros(1,nmodeE);     %thetaֵ
heffE=zeros(1,nmodeE);      %heffֵ
AE=zeros(1,nmodeE);         %Aֵ
phiE=zeros(1,nmodeE);       %���Ŀռ�ֲ����ʽ��cos(phi)�е�phiֵ
sumPsubE=zeros(1,nmodeE);   %�ĵײ���ܹ���
sumPcoreE=zeros(1,nmodeE);  %о����ܹ���
sumPcoverE=zeros(1,nmodeE); %���ǲ���ܹ���

%���TEģ�����������Ŀռ�ֲ��������Ŀռ�ֲ�
%����������Ĳ������ǵ�����PPT���壨ǰ�벿�����ݣ��ı��ʽ�ĸ�������
%�������еĳ��Ŀռ�ֲ����ǵ�����ѧ�α��еı��ʽ
%���ռ�ֲ��ֶΣ�-inf,0��,(0,d),(d,inf)
for j=1:nmodeE
    hE(j)=solutionE(1,j)/d;
    pE(j)=sqrt((n1^2-n2^2)*k^2-(solutionE(1,j)/d)^2);
    qE(j)=sqrt((n1^2-n3^2)*k^2-(solutionE(1,j)/d)^2);
    betaE(j)=sqrt(n1^2*k^2-(solutionE(1,j)/d)^2);
    thetaE(j)=asind(betaE(j)/(k*n1));
    heffE(j)=d+(1/pE(j))+(1/qE(j));
    AE(j)=2*hE(j)*sqrt(w0*u0/(heffE(j)*betaE(j)*(hE(j)^2+pE(j)^2)));
    phiE(j)=atan(pE(j)/hE(j));
    Esub=AE(j)*exp(pE(j)*X);                                            %TEģ�ĳĵײ�糡�ռ�ֲ�
    Ecore=AE(j)*cos(hE(j)*X-phiE(j))/cos(phiE(j));                      %TEģ��о��糡�ռ�ֲ�
    Ecover=(AE(j)*cos(hE(j)*d-phiE(j))/cos(phiE(j)))*exp(-qE(j)*(X-d)); %TEģ�ĸ��ǲ�糡�ռ�ֲ�
    PsubE=betaE(j)*(Esub^2)/(2*w0*u0);                                  %TEģ�ĳĵײ�糡�ռ�ֲ�
    PcoreE=betaE(j)*(Ecore^2)/(2*w0*u0);                                %TEģ��о��糡�ռ�ֲ�
    PcoverE=betaE(j)*(Ecover^2)/(2*w0*u0);                              %TEģ�ĸ��ǲ�糡�ռ�ֲ�
    
    %��ʾTEģ�ĸ���糡�������Ŀռ�ֲ�
    fprintf('��%d��TE��ģ\n',j);
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

%���TEģ������ܹ���
for j=1:nmodeE
    sumPsubE(j)=AE(j)^2*betaE(j)/(4*w0*u0*pE(j));
    sumPcoreE(j)=AE(j)^2*betaE(j)*(hE(j)^2+pE(j)^2)*(d+pE(j)/(hE(j)^2+pE(j)^2)+qE(j)/(hE(j)^2+qE(j)^2))/(4*w0*u0*hE(j)^2);
    sumPcoverE(j)=AE(j)^2*betaE(j)*(hE(j)^2+pE(j)^2)/(4*w0*u0*qE(j)*(hE(j)^2+qE(j)^2));
end

%��ʾTEģ������
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

%��ʾTEģ�ĵ糡�������Ŀռ�ֲ�
X1=-d:0.001*d:2*d;
for j=1:nmodeE
    Esub=AE(j)*exp(pE(j)*X1);
    Ecore=AE(j)*cos(hE(j)*X1-phiE(j))/cos(phiE(j));
    Ecover=(AE(j)*cos(hE(j)*d-phiE(j))/cos(phiE(j)))*exp(-qE(j)*(X1-d));
    PsubE=betaE(j)*power(Esub,2)/(2*w0*u0);
    PcoreE=betaE(j)*power(Ecore,2)/(2*w0*u0);
    PcoverE=betaE(j)*power(Ecover,2)/(2*w0*u0);
    PicE=Esub.*(X1<0)+Ecore.*((0<X1)&(X1<d))+Ecover.*(X1>d);        %�糡�ķֶκ���
    PicPE=PsubE.*(X1<0)+PcoreE.*((0<X1)&(X1<d))+PcoverE.*(X1>d);    %�������ķֶκ���
    
    %�������������ͼ�ηָ��ʾ
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



%*******************TM����ģ�������*****************
%*******************ע����TEģ����********************
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

%��������
for j=1:nmodeM
    hM(j)=solutionM(1,j)/d;
    pM(j)=sqrt((n1^2-n2^2)*k^2-(solutionM(1,j)/d)^2);
    qM(j)=sqrt((n1^2-n3^2)*k^2-(solutionM(1,j)/d)^2);
    betaM(j)=sqrt(n1^2*k^2-(solutionM(1,j)/d)^2);
    thetaM(j)=asind(betaM(j)/(k*n1));
    heffM(j)=d+(n1^2*n3^2)*(hM(j)^2+qM(j)^2)/((n3^4*hM(j)^2+n1^4*qM(j)^2)*qM(j))+(n1^2*n2^2)*(hM(j)^2+pM(j)^2)/((n2^4*hM(j)^2+n1^4*pM(j)^2)*pM(j));
    AM(j)=2*hM(j)*sqrt((n1^2*n2^4)*w0*eps0/(betaM(j)*heffM(j)*(n2^4*hM(j)^2+n1^4*pM(j)^2)));
    phiM(j)=atan(n1^2*pM(j)/(n2^2*hM(j)));
    Hsub=AM(j)*exp(pM(j)*X);                                            %�ĵײ�ų�����
    Hcore=AM(j)*cos(hM(j)*X-phiM(j))/cos(phiM(j));                      %о��ų�����
    Hcover=(AM(j)*cos(hM(j)*d-phiM(j))/cos(phiM(j)))*exp(-qM(j)*(X-d)); %���ǲ�ų�����
    PsubM=betaM(j)*(Hsub^2)/(2*w0*eps0*n2^2);                           %�ĵײ�����������
    PcoreM=betaM(j)*(Hcore^2)/(2*w0*eps0*n1^2);                         %о������������
    PcoverM=betaM(j)*(Hcover^2)/(2*w0*eps0*n3^2);                       %���ǲ�����������
    sumPsubM(j)=vpa(int(PsubM,X,-inf,0),5);                             %�ĵײ�������
    sumPcoreM(j)=vpa(int(PcoreM,X,0,d),5);                              %о��������
    sumPcoverM(j)=vpa(int(PcoverM,X,d,inf),5);                          %���ǲ�������
    
    %��ʾ�ų����������Ŀռ�ֲ�����
    fprintf('��%d��TM��ģ\n',j);
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

%��ʾ�ų����������Ŀռ�ֲ�ͼ��
X2=-d:0.001*d:2*d;
for j=1:3
    Hsub=AM(j)*exp(pM(j)*X2);
    Hcore=AM(j)*cos(hM(j)*X2-phiM(j))/cos(phiM(j));
    Hcover=(AM(j)*cos(hM(j)*d-phiM(j))/cos(phiM(j)))*exp(-qM(j)*(X2-d));
    PsubM=betaM(j)*power(Hsub,2)/(2*w0*eps0*n2^2);
    PcoreM=betaM(j)*power(Hcore,2)/(2*w0*eps0*n1^2);
    PcoverM=betaM(j)*power(Hcover,2)/(2*w0*eps0*n3^2);
    PicM=Hsub.*(X2<0)+Hcore.*(X2>0&X2<d)+Hcover.*(X2>d);        %�ų��ֲ�����
    PicPM=PsubM.*(X2<0)+PcoreM.*(X2>0&X2<d)+PcoverM.*(X2>d);    %�������ֲ�����
    
    %�������������ͼ�ηָ��ʾ
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

    
%*********��ʾTMģ������**********

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


%**************************���ݵ�������*************************
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
%******************************���ݵ������**********************
