clear all;clc;
tic
global N;
global ecm;
global ecm1;
global scm;
global scm1;
global D;
global D1;
global G;
global Fr0;
%global gradE;
% N particles
N=1;
ecm=-0.1092;
ecm1=-0.067;
scm=-0.5;
scm1=-0.3333;
D=5.152;
D1=5.6054;
G=4;
Fr0=-1;
%gradE=[1,0,0;0,0,0;0,0,-1];


%initial conditions
%r=r0
r=[2;0;6];
Q=zeros(3,3);
Om=zeros(3,1);
sp=1e-3;
P=sp*(1-2*rand(3,1));


%data and time
dt=0.005;
tmax=400;
time=0:dt:tmax-dt;
tmax0=tmax/dt;
datar=zeros(3,tmax0);%location
dataP=zeros(3,tmax0);%dipole moment
dataQ=zeros(6,tmax0);%quadrupole moment
dataOm=zeros(3,tmax0);%rotational rate
datar0=zeros(1,tmax0);%distance from 0,0


%save initial data
datar(:,1)=r(:);
dataP(:,1)=P(:);
dataQ(:,1)=[Q(1,1);Q(1,2);Q(1,3);Q(2,2);Q(2,3);Q(3,3)];
dataOm(:,1)=Om(:);
datar0(1)=sqrt(r(1)^2+r(2)^2+r(3)^2);

%loop
for t=2:tmax0
    [dp,dq] = fpq1(P,Q,r,Om);
    P = P + dt*dp;
    Q = Q + dt*dq;
    dr = fr1(P,Q,r);
    r = r + dt*dr;
    Om = fom1(P,Q,r);
    
    
    % Save data
    datar(:,t)=r(:);
    dataP(:,t)=P(:);
    dataQ(:,t)=[Q(1,1);Q(1,2);Q(1,3);Q(2,2);Q(2,3);Q(3,3)];
    dataOm(:,t)=Om(:);
    datar0(t)=sqrt(r(1)^2+r(2)^2+r(3)^2);
end



% Plot
figure(1);
plot3(datar(1,:),datar(2,:),datar(3,:));grid on;hold on;axis equal;
plot3(datar(1,(tmax-200)/dt:tmax0),datar(2,(tmax-200)/dt:tmax0),datar(3,(tmax-200)/dt:tmax0),'LineWidth',1);grid on;hold on;axis equal;
scatter3(datar(1,1),datar(2,1),datar(3,1),'*');
scatter3(datar(1,tmax0),datar(2,tmax0),datar(3,tmax0),'filled');
xlabel('x');ylabel('y');zlabel('z');
% plot(datar(1,:),datar(3,:));grid on;hold on;axis equal;
% plot(datar(1,(tmax-150)/dt:tmax0),datar(3,(tmax-150)/dt:tmax0),'LineWidth',1);grid on;hold on;axis equal;
% scatter(datar(1,1),datar(3,1),'*');
% scatter(datar(1,tmax0),datar(3,tmax0),'filled');
% xlabel('x');ylabel('z');

figure(2);
subplot(1,2,1);
plot(time,dataOm(2,:));grid on;
xlabel('t');ylabel('\Omega');

subplot(1,2,2);
plot(time,datar0);grid on;
xlabel('t');ylabel('Distance from (0,0)');

toc