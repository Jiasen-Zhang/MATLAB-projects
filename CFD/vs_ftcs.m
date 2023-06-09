clear all;clc;
tic
%参数设置
err = 1e-4;
re=100;
L=1;               %空间几何尺寸
n=50;
dh=L/n;
dt=1e-3;
psi=zeros(n,n);%流函数
w=zeros(n,n);
w0=w;%
u=zeros(n,n);
v=u;

%
I = speye(n-2);
e = ones(n-2,1);
T = spdiags([e -4*e e],[-1 0 1],n-2,n-2);
S = spdiags([e e],[-1 1],n-2,n-2);
A = (kron(I,T) + kron(S,I)) / dh^2;

%
for j=1:n
   v(n,j)=j*dh*(1-j*dh);
end

for k=1:1000000
    
    % boundary
    w(:,1)=-2*(psi(:,2)-psi(:,1))/(dh^2);
    w(:,n)=-2*(psi(:,n-1)-psi(:,n))/(dh^2);
    for j=1:n
        w(n,j)=-2*(psi(n-1,j)-psi(n,j))/(dh^2) + j*dh*(1-j*dh)*2/dh;
    end
    w(1,:)=-2*(psi(2,:)-psi(1,:))/(dh^2);
    
    % corner
    w(1,1)=(w(1,2)+w(2,1))/2;
    w(1,n)=(w(1,n-1)+w(2,n))/2;
    w(n,1)=(w(n,2)+w(n-1,1))/2;
    w(n,n)=(w(n,n-1)+w(n-1,n))/2;
    
    % FTCS
    u(2:n-1,2:n-1)=(psi(2:n-1,3:n)-psi(2:n-1,1:n-2))/(2*dh);
    v(2:n-1,2:n-1)=-(psi(3:n,2:n-1)-psi(1:n-2,2:n-1))/(2*dh);
    uw = u.*w;
    vw = v.*w;

    % FTCS
    dw = -(uw(3:n,2:n-1)-uw(1:n-2,2:n-1))/(2*dh)...
        -(vw(2:n-1,3:n)-vw(2:n-1,1:n-2))/(2*dh)...
        +(w(3:n,2:n-1)+w(1:n-2,2:n-1)+w(2:n-1,3:n)+w(2:n-1,1:n-2)-4*w(2:n-1,2:n-1))/(re*dh^2);
%     dw2 = -u(2:n-1,2:n-1).*(w(3:n,2:n-1)-w(1:n-2,2:n-1))/(2*dh)...
%          -v(2:n-1,2:n-1).*(w(2:n-1,3:n)-w(2:n-1,1:n-2))/(2*dh)...
%          +(w(3:n,2:n-1)+w(1:n-2,2:n-1)+w(2:n-1,3:n)+w(2:n-1,1:n-2)-4*w(2:n-1,2:n-1))/(re*dh^2);

    dw = dt*dw;
    w(2:n-1,2:n-1) = w(2:n-1,2:n-1) + dw;
    diff = max(max(abs(dw)));

    
    % psi: 0 at boundaries
    wvec = reshape(w(2:n-1,2:n-1),(n-2)^2,1);
    psiv = -A\wvec;
    psi(2:n-1,2:n-1) = reshape(psiv,n-2,n-2);
    
    
    if mod(k,500)==0
        fprintf('%d\t %e\n',k,diff);
    end
    if diff<err
        break;
    end
end

figure(1);
[x,y] = meshgrid(0:dh:L-dh,0:dh:L-dh);axis equal;
streamslice(x,y,v,u,10);
figure(2);
contour(x,y,psi,'ShowText','on');axis equal;
figure(3);
subplot(1,2,1);
contour(x,y,u,[-0.07:0.01:0.07],'ShowText','on');axis equal;
subplot(1,2,2);
contour(x,y,v,[-0.04:0.02:0.22],'ShowText','on');axis equal;
toc