clear all;clc;
tic
%参数设置
%D为矩阵
e=1e-5;
re=3000;
L=1;               %空间几何尺寸
n=50;
dh=L/n;            %dx=dy=dh
dt=1/500;           %时间步长

u=zeros(n,n);
v=u;p=u;u0=u;v0=v;
D=u;

%初始条件
p(:,:)=1;
%边界条件
for j=1:n
    v(1,j)=0.5;%dh*j*(1-dh*j);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:10000
    
    %求散度
    fux=u(2:n-1,2:n-1).*(u(2:n-1,2:n-1)-u(1:n-2,2:n-1))/dh;
    fuy=v(2:n-1,2:n-1).*(u(2:n-1,2:n-1)-u(2:n-1,1:n-2))/dh;
    fvx=u(2:n-1,2:n-1).*(v(2:n-1,2:n-1)-v(1:n-2,2:n-1))/dh;
    fvy=v(2:n-1,2:n-1).*(v(2:n-1,2:n-1)-v(2:n-1,1:n-2))/dh;
    visx=(u(3:n,2:n-1)+u(1:n-2,2:n-1)+u(2:n-1,3:n)+u(2:n-1,1:n-2)-4*u(2:n-1,2:n-1))/(re*dh*dh);
    visy=(v(3:n,2:n-1)+v(1:n-2,2:n-1)+v(2:n-1,3:n)+v(2:n-1,1:n-2)-4*v(2:n-1,2:n-1))/(re*dh*dh);
    
    u0(2:n-1,2:n-1) = ((p(2:n-1,2:n-1)-p(3:n,2:n-1))/dh+visx-fux-fuy)*dt+u(2:n-1,2:n-1);
    v0(2:n-1,2:n-1) = ((p(2:n-1,2:n-1)-p(2:n-1,3:n))/dh+visy-fvx-fvy)*dt+v(2:n-1,2:n-1); 
    D(2:n-1,2:n-1)=(u0(2:n-1,2:n-1)-u0(1:n-2,2:n-1)+v0(2:n-1,2:n-1)-v0(2:n-1,1:n-2))/dh;

    

    %检验散度是否为零
    while(max(max((abs(D))))>1e-4)   
        dp2 = -D(2:n,2:n)*dh/4;
        dp = dp2*dh/dt;
        u0(2:n,2:n)=u0(2:n,2:n)+dp2;
        u0(1:n-1,2:n)=u0(1:n-1,2:n)-dp2;
        v0(2:n,2:n)=v0(2:n,2:n)+dp2;
        v0(2:n,1:n-1)=v0(2:n,1:n-1)-dp2;
        D(2:n,2:n)=(u0(2:n,2:n)-u0(1:n-1,2:n)+v0(2:n,2:n)-v0(2:n,1:n-1))/dh;
        p(2:n,2:n)=p(2:n,2:n)+dp; 
    end
%     while(max(max((abs(D))))>1e-4)        
%         for j=2:n
%             for i=2:n
%                 dp=-D(i,j)*dh*dh/(4*dt);
%                 dp2=dt*dp/dh;
%                 u0(i,j)=u0(i,j)+dp2;
%                 u0(i-1,j)=u0(i-1,j)-dp2;
%                 v0(i,j)=v0(i,j)+dp2;
%                 v0(i,j-1)=v0(i,j-1)-dp2;
%                 D(i,j)=(u0(i,j)-u0(i-1,j)+v0(i,j)-v0(i,j-1))/dh;
%                 p(i,j)=p(i,j)+dp;             
%             end   
%         end
%     end
    
    % boundary condition
    u0(:,1)=0;
    u0(:,n)=0;
    v0(n,:)=0;
    for j=1:n
        v0(1,j)=0.5;%dh*j*(1-dh*j);
    end
    
        
    %检验收敛
    du=max(max(abs(u0-u)));
    dv=max(max(abs(v0-v)));
    u=u0;v=v0;
    err=max(du,dv);
    if mod(k,500)==0
        fprintf('%d\t %e\n',k,err);
    end
    if err<=e
        break;
    end
    
end
%计算完毕
figure(1);
streamslice(v,u,10);

%计算流函数
psi=zeros(n,n);
w=(v(3:n,2:n-1)-v(1:n-2,2:n-1))/(2*dh)-(u(2:n-1,3:n)-u(2:n-1,1:n-2))/(2*dh);
for k=1:100000
    err1=(psi(3:n,2:n-1)+psi(1:n-2,2:n-1)+psi(2:n-1,3:n)+psi(2:n-1,1:n-2)+w*dh*dh)/4-psi(2:n-1,2:n-1);
    psi(2:n-1,2:n-1)=psi(2:n-1,2:n-1)+err1;  
    if max(max(abs(err1)))<1e-8
        break;
    end
end
figure(2);
contour(psi,n+1);
toc