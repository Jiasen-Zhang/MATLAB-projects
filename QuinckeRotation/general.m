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
N=1;
ecm=-0.1092;
ecm1=-0.067;
scm=-0.5;
scm1=-0.3333;
D=5.152;
D1=5.6054;
G=3
Fr0=-1;



%initial conditions
r0=[2,2,2,2,-2,-2,-2,-2;                                                                                                                                                                         
    0,2,-2,-2,2,2,-2,-2;
    2,-2,2,-2,2,-2,2,-2];
Q=zeros(3*N,3);
r=zeros(3*N,1);
Om=zeros(3*N,1);
sp=1e-3;
P=sp*(1-2*rand(3*N,1));
for i=1:N
    xi=3*i-2;
    zi=3*i;
    r(xi:zi)=r0(:,i);
    %P(xi+1)=0;
end


%data and time
dt=0.005;
tmax=700;
time=0:dt:tmax-dt;
tmax0=tmax/dt;
datar=zeros(3*N,tmax0);%location
dataP=zeros(3*N,tmax0);%dipole moment
dataQ=zeros(6*N,tmax0);%quadrupole moment
dataOm=zeros(3*N,tmax0);%rotational rate
dataom3=zeros(N,tmax0);
datar0=zeros(N,tmax0);%distance between i and (0,0)


%save initial data
for i=1:N
    xi=3*i-2;
    yi=3*i-1;
    zi=3*i;
    datar(xi:zi,1)=r(xi:zi);
    dataP(xi:zi,1)=P(xi:zi);
    dataQ(6*i-5:6*i,1)=[Q(xi,1);Q(xi,2);Q(xi,3);Q(yi,2);Q(yi,3);Q(zi,3)];
    dataOm(xi:zi,1)=Om(xi:zi);
    dataom3(i,1)=sqrt( Om(xi)^2+Om(yi)^2+Om(zi)^2 );
    datar0(i,1)=sqrt( r(xi)^2+r(yi)^2+r(zi)^2 );
end



%loop
for t=2:tmax0
    [dp,dq] = fpqn(P,Q,r,Om);
    P = P + dt*dp;
    Q = Q + dt*dq;
    dr = frn(P,Q,r);
    r = r + dt*dr;
    Om = fomn(P,Q,r);
    % Save data
    for i=1:N
        xi=3*i-2;
        yi=3*i-1;
        zi=3*i;
        datar(xi:zi,t)=r(xi:zi);
        dataP(xi:zi,t)=P(xi:zi);
        dataQ(6*i-5:6*i,t)=[Q(xi,1);Q(xi,2);Q(xi,3);Q(yi,2);Q(yi,3);Q(zi,3)];
        dataOm(xi:zi,t)=Om(xi:zi);
        dataom3(i,t)=sqrt( Om(xi)^2+Om(yi)^2+Om(zi)^2 );
        datar0(i,t)=sqrt( r(xi)^2+r(yi)^2+r(zi)^2 );

    end

end



% Plot
t0=100;
for i=1:N
     xi=3*i-2;
     yi=3*i-1;
     zi=3*i;
     subplot(1,3,1);
      grid on;hold on;axis equal;
%      plot3(datar(xi,:),datar(yi,:),datar(zi,:));
%      plot3(datar(1,(tmax-t0)/dt:tmax0),datar(2,(tmax-t0)/dt:tmax0),datar(3,(tmax-t0)/dt:tmax0),'LineWidth',1);
%       scatter3(datar(xi,1),datar(yi,1),datar(zi,1),'*');
%       scatter3(datar(xi,tmax0),datar(yi,tmax0),datar(zi,tmax0),'filled');
%       xlabel('x');ylabel('y');zlabel('z');
     % 
     plot(datar(xi,:),datar(zi,:),'b');
     plot(datar(xi,(tmax-t0)/dt:tmax0),datar(zi,(tmax-t0)/dt:tmax0),'r','Linewidth',2);
     scatter(datar(xi,1),datar(zi,1),'*','r');
     %scatter(datar(xi,tmax0),datar(zi,tmax0),'filled','b');
     xlabel('x');ylabel('y');
     title('(a)');

     subplot(1,3,3);
     title('ROTATION');
     plot(time,dataom3(i,:));grid on;hold on;
     xlabel('t');ylabel('\Omega');
%      
%      subplot(2,2,3);
%      plot(time,datar0(i,:));grid on;hold on;
%      xlabel('t');ylabel('distance');
     
     subplot(1,3,2);axis equal;grid on;hold on;
     %scatter(datar(xi,tmax0),datar(zi,tmax0),'filled','b');
     plot(datar(xi,(tmax-t0)/dt:tmax0),datar(zi,(tmax-t0)/dt:tmax0),'r');
     xlabel('x');ylabel('y');
     title('(b)');
end

toc