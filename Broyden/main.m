function main
%13307110108 �ż�ɭ
%�öγ������ڼ��㲻ͬphi�µ����Ž⣬��Ӧ����ֵ�͵�������
%��������ú�����Сֵ�����phi��ɢ��ͼ
clc;
i=1;
f=@(x)0.5*x(1)^2+x(2)^2-x(1)*x(2);%���ú������ʽ
x0=[2;2];
e=1e-05;
%������������
best=eye(2,10);
fmini=eye(1,10);
ki=eye(1,10);
for phi=0:0.1:1 %phi��0��1�����Ϊ0.1
    [ bestsolution,fmin,k ] =broyden(f,phi,x0,e);
    best(:,i)=bestsolution;
    fmini(i)=fmin;
    ki(i)=k;
    scatter(phi,fmin,10,'filled');%������Сֵfmin��phi�Ĺ�ϵɢ��ͼ
    hold on;
    i=i+1;
end
best%�����г���phiֵ�µ����Ž�
fmini%�����г���phiֵ�µ�fmin
ki%�����г���phiֵ�µĵ�������
end

