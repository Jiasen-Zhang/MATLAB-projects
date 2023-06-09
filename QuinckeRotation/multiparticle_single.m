clear all;clc;
tic
% N particles
N=1;
%parameters
ecm=-0.1092;
ecm1=-0.067;
scm=-0.5;
scm1=-0.3333;
D=5.152;
D1=5.6054;
G=0.4;
Fr0=0;
gradE=[1,0,0;0,0,0;0,0,-1];
%
err=1;
dt=0.01;
%initial conditions
r0=[2.5,11,4.5;0,0,0;6,19,6];
r=zeros(3*N,1);
P=zeros(3*N,1);
P(1:3)=-((1-2*rand(3,1))*4.5e-4+5.5e-4);% from -1e-4 to -1e-3
P(4:6)=-((1-2*rand(3,1))*4.5e-4+5.5e-4);% from -1e-4 to -1e-3

Q=zeros(3*N,3);
Om=zeros(3*N,1);
for i=1:N
    xi=3*i-2;
    yi=3*i-1;
    zi=3*i;
    r(xi:zi)=r0(:,i);
    P(yi)=0;
end
%
Fel=zeros(3*N,1);
Tel=zeros(3*N,1);
%
tmax=200;
tmax0=tmax/dt;
datar=zeros(3*N,tmax);
dataP=zeros(3*N,tmax);
dataQ=zeros(6*N,tmax);
dataOm=zeros(3*N,tmax);
%
for i=1:N
    xi=3*i-2;
    yi=3*i-1;
    zi=3*i;
    datar(xi:zi,1)=r(xi:zi);
    dataP(xi:zi,1)=P(xi:zi);
    dataQ(6*i-5:6*i,1)=[Q(xi,1);Q(xi,2);Q(xi,3);Q(yi,2);Q(yi,3);Q(zi,3)];
    dataOm(xi:zi,1)=Om(xi:zi);
end
% %
for t=2:tmax0
    % Get Rij, R0ij, Epij, Eqij, Gpij
    % Get Pi, Qi
    for i=1:N
        xi=3*i-2;
        zi=3*i;
        Ei=[r(xi);0;-r(zi)];
        %
        dPi1=P(xi:zi)+G*ecm*(-Ei);
        dPi2=(P(xi:zi)+G*scm*(-Ei))/D;
        dPi=cross(Om(xi:zi),dPi1)-dPi2;
        P(xi:zi)=P(xi:zi)+dt*dPi;
        %
        dQi1=Q(xi:zi,:)+2*G*ecm1*(-gradE);
        dQi2=[cross(Om(xi:zi),dQi1(1,:));cross(Om(xi:zi),dQi1(2,:));cross(Om(xi:zi),dQi1(3,:))];
        dQi3=(Q(xi:zi,:)+2*G*scm1*(-gradE))/D1;
        dQi=dQi2+dQi2'-dQi3;
        Q(xi:zi,:)=Q(xi:zi,:)+dt*dQi;
    end
    % Get Fd2i0, Frep, Td4i0
    % Get Fel and Tel
    for i=1:N
        xi=3*i-2;   
        yi=3*i-1;
        zi=3*i;
        Ei=[r(xi);0;-r(zi)];
        %
        Fd1i=-4*pi*G*(-gradE)*P(xi:zi);
        Fd2i=0;
        Fd3i=0;
        Fel(xi:zi)=Fd1i+Fd2i+Fd3i;
        %
        Td1i=-4*pi*G*(cross(P(xi:zi),-Ei)+[Q(yi,3);-Q(zi,1)-Q(xi,3);Q(yi,1)]);
        Td2i=0;
        Td3i=0;
        Td4i=0;
        Tel(xi:zi)=Td1i+Td2i+Td3i+Td4i;
    end
    % Get dui1, dui2, dOmi1
    % Get ri, Omi
    for i=1:N
        xi=3*i-2;
        zi=3*i;
        dui=Fel(xi:zi)/(6*pi);
        Om(xi:zi)=Tel(xi:zi)/(8*pi);
        r(xi:zi)=r(xi:zi)+dt*dui;
    end
    % Save data
    for i=1:N
        xi=3*i-2;
        yi=3*i-1;
        zi=3*i;
        datar(xi:zi,t)=r(xi:zi);
        dataP(xi:zi,t)=P(xi:zi);
        dataQ(6*i-5:6*i,t)=[Q(xi,1);Q(xi,2);Q(xi,3);Q(yi,2);Q(yi,3);Q(zi,3)];
        dataOm(xi:zi,t)=Om(xi:zi);
    end
    
end
% Plot
for i=1:N
     xi=3*i-2;
     yi=3*i-1;
     zi=3*i;
    plot(datar(xi,:),datar(zi,:));grid on;hold on;
    xlabel('x');ylabel('z');
end
axis equal;
toc