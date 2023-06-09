
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
% N particles
N=2;
ecm=-0.1092;
ecm1=-0.067;
scm=-0.5;
scm1=-0.3333;
D=5.152;
D1=5.6054;
G=1.3;
Fr0=-1;


%initial conditions
r0=[2,-2,2,-2;
    0,0,0,2;
    6,-2,-2,-2];
r=zeros(3*N,1);
Q=zeros(3*N,3);
Om=zeros(3*N,1);
sp=1e-3;
P=sp*(1-2*rand(3*N,1));
for i=1:N
    xi=3*i-2;
    zi=3*i;
    r(xi:zi)=r0(:,i);
    P(xi+1)=0;
end
%
%
tmax=400;

%  assemble matrix
data0=zeros(18*N,1); %r,P and Q
data0(1:3*N,1)=r;
data0(3*N+1:6*N,1)=P;
data0(6*N+1:9*N)=Q(:,1);
data0(9*N+1:12*N)=Q(:,2);
data0(12*N+1:15*N)=Q(:,3);

[t,output]=ode15s(@F,[0,tmax],data0);

% get the data
datar = output(:,1:3*N)';
dataom = output(:,15*N+1:18*N)';
datap = output(:,3*N+1:6*N)';

%plot
tmax0 = length(t);
for i=1:N
     xi=3*i-2;
     yi=3*i-1;
     zi=3*i;
     subplot(2,2,1);
%      plot3(datar(xi,:),datar(yi,:),datar(zi,:));grid on;hold on;axis equal;
%      scatter3(datar(xi,1),datar(yi,1),datar(zi,1),'*');
%      scatter3(datar(xi,tmax0),datar(yi,tmax0),datar(zi,tmax0),'filled');
%      xlabel('x');ylabel('y');zlabel('z');
     plot(datar(xi,:),datar(zi,:));grid on;hold on;axis equal;
     plot(datar(xi,tmax0-200:tmax0),datar(zi,tmax0-200:tmax0));
     scatter(datar(xi,1),datar(zi,1),'*');
     scatter(datar(xi,tmax0),datar(zi,tmax0),'filled');
     xlabel('x');ylabel('z');

     subplot(2,2,2);
     plot(t,dataom(yi,:));grid on;hold on;
     xlabel('t');ylabel('\Omega');
end

toc
