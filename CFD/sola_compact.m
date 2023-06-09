clear all;clc;
tic
%参数设置
%D为矩阵
e=1e-4;
re=100;
L=1;               %空间几何尺寸
n=50;
dh=L/n;            %dx=dy=dh
dt=1e-3;           %时间步长

u=zeros(n,n);
v=u;p=u;u0=u;v0=v;
D=u;
erru=D;errv=D;
err=D;

%初始条件
p(:,:)=1;

%边界条件
for j=1:n
    v(1,j)=dh*j*(1-dh*j);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:100000
    
    vj1=compactd(v,1,1,dh);  %v j方向1阶导数
    vj2=compactd(v,1,2,dh);  %v j方向2阶导数
    vi1=compactd(v,2,1,dh);  %v i方向1阶导数
    vi2=compactd(v,2,2,dh);  %v i方向2阶导数 
    uj1=compactd(u,1,1,dh);  %u j方向1阶导数    
    uj2=compactd(u,1,2,dh);  %u j方向2阶导数    
    ui1=compactd(u,2,1,dh);  %u i方向1阶导数 
    ui2=compactd(u,2,2,dh);  %u i方向2阶导数
    pj1=compactd(p,1,1,dh);  %p j方向1阶导数
    pi1=compactd(p,2,1,dh);  %p i方向1阶导数
    fui = u.*ui1;
    fuj = v.*uj1;
    visu = (ui2+uj2)/re;
    fvi = u.*vi1;
    fvj = v.*vj1;
    visv = (vi2+vj2)/re;
    
    
    err(2:n-1,2:n-1) = -pi1(2:n-1,2:n-1) ...
        +visu(2:n-1,2:n-1)-fui(2:n-1,2:n-1)-fuj(2:n-1,2:n-1);
    u0 = err*dt + u;
    err(2:n-1,2:n-1) = -pj1(2:n-1,2:n-1) ...
        +visv(2:n-1,2:n-1)-fvi(2:n-1,2:n-1)-fvj(2:n-1,2:n-1);         
    v0 = err*dt +v;
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
    end
    
    %边界条件
    % boundary condition
    u0(:,1)=0;v0(:,1)=0;
    u0(:,n)=0;v0(:,n)=0;
    v0(n,:)=0;u0(n,:)=0;
    u0(1,:)=0;
    for j=1:n
        v0(1,j)=dh*j*(1-dh*j);
    end
    
    %检验收敛
    du=max(max(abs(u0-u)));
    dv=max(max(abs(v0-v)));
    errk=max(du,dv);
    u=u0;v=v0;
    if mod(k,1)==0
        fprintf('%d\t %e\n',k,errk);
    end
    if errk<=e
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
    psi0 = (psi(3:n,2:n-1)+psi(1:n-2,2:n-1)+psi(2:n-1,3:n)+psi(2:n-1,1:n-2)+w*(dh^2))/4;
    err1=psi0-psi(2:n-1,2:n-1);
    psi(2:n-1,2:n-1)=psi0;  
    if max(max(abs(err1)))<1e-8
        break;
    end
end
figure(2);
contour(psi,n);
toc

function d = compactd(u,z,m,dh)
%   4阶紧致差分
%   u求导的矩阵
%   z为方向i(2)或j(1)
%   m为求导阶数1或2
[ni,nj]=size(u);
dd=zeros(ni,nj);

if z==2
    %i方向一阶导数
    if m==1
        Fi=zeros(ni,ni);
        Fi(1,1)=1;Fi(ni,ni)=1;
        Fi(1,2)=3;Fi(ni,ni-1)=3;
        for i=2:ni-1
            Fi(i,i)=2/3;
            Fi(i,i-1)=1/6;
            Fi(i,i+1)=1/6;
        end
        dd(1,:)=(-17*u(1,:)+9*u(2,:)+9*u(3,:)-u(4,:))/(6*dh);
        dd(ni,:)=(u(ni-3,:)-9*u(ni-2,:)-9*u(ni-1,:)+17*u(ni,:))/(6*dh);
        dd(2:ni-1,:)=(u(3:ni,:)-u(1:ni-2,:))/(2*dh);
        d=Fi\dd;
    end
    
    %i方向二阶导数
    if m==2
        Si=zeros(ni,ni);
        Si(1,1)=1;Si(ni,nj)=1;
        Si(1,2)=11;Si(ni,ni-1)=11;
        for i=2:ni-1
            Si(i,i)=5/6;
            Si(i,i-1)=1/12;
            Si(i,i+1)=1/12;
        end
        dd(1,:)=(13*u(1,:)-27*u(2,:)+15*u(3,:)-u(4,:))/(dh*dh);
        dd(ni,:)=(-u(ni-3,:)+15*u(ni-2,:)-27*u(ni-1,:)+13*u(ni,:))/(dh*dh);
        dd(2:ni-1,:)=(u(3:ni,:)-2*u(2:ni-1,:)+u(1:ni-2,:))/(dh*dh);
        d=Si\dd;
    end
end

if z==1
    %j方向一阶导数
    if m==1
        Fj=zeros(nj,nj);
        Fj(1,1)=1;Fj(nj,nj)=1;
        Fj(2,1)=3;Fj(nj-1,nj)=3;
        for i=2:nj-1
            Fj(i,i)=2/3;
            Fj(i-1,i)=1/6;
            Fj(i+1,i)=1/6;
        end
        dd(:,1)=(-17*u(:,1)+9*u(:,2)+9*u(:,3)-u(:,4))/(6*dh);
        dd(:,nj)=(u(:,nj-3)-9*u(:,nj-2)-9*u(:,nj-1)+17*u(:,nj))/(6*dh);
        dd(:,2:nj-1)=(u(:,3:nj)-u(:,1:nj-2))/(2*dh);
        d=dd/Fj;
    end
    %j方向二阶导数
    if m==2
        Sj=zeros(nj,nj);
        Sj(1,1)=1;Sj(nj,nj)=1;
        Sj(2,1)=11;Sj(nj-1,nj)=11;
        for i=2:nj-1
            Sj(i,i)=5/6;
            Sj(i-1,i)=1/12;
            Sj(i+1,i)=1/12;
        end
        dd(:,1)=(13*u(:,1)-27*u(:,2)+15*u(:,3)-u(:,4))/(dh*dh);
        dd(:,nj)=(-u(:,nj-3)+15*u(:,nj-2)-27*u(:,nj-1)+13*u(:,nj))/(dh*dh);
        dd(:,2:nj-1)=(u(:,3:nj)-2*u(:,2:nj-1)+u(:,1:nj-2))/(dh*dh);
        d=dd/Sj;
    end
end

end
