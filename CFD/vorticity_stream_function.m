clear all;clc;
tic
%参数设置
re=1000;
L=1;               %空间几何尺寸
n=50;
dh=L/n;
dt=1e-3;           %时间步长
psi=zeros(n+1,n+1);%流函数
xi=zeros(n+1,n+1); %涡量
u=zeros(n,n);v=u;
rho=1;             %密度?
%
for k=1:100000
    err=0;
    %边界条件
    for i=2:n
        xi(i,1)=-2*(psi(i,2)-psi(i,1))/(dh^2);
        xi(i,n+1)=-2*(psi(i,n)-psi(i,n+1))/(dh^2);
    end
    for j=2:n
        xi(1,j)=-2*(psi(2,j)-psi(1,j)+dh)/(dh^2);
        xi(n+1,j)=-2*(psi(n,j)-psi(n+1,j))/(dh^2);
    end
    xi(1,1)=(xi(1,2)+xi(2,1))/2;
    xi(1,n+1)=(xi(1,n)+xi(2,n+1))/2;
    xi(n+1,1)=(xi(n+1,2)+xi(n,1))/2;
    xi(n+1,n+1)=(xi(n+1,n)+xi(n,n+1))/2;
    %控制方程
    for i=2:n
        for j=2:n
            %u v
            u(i,j)=(psi(i,j+1)-psi(i,j-1))/(2*dh);
            v(i,j)=-((psi(i+1,j)-psi(i-1,j))/(2*dh));
            %stream function
            err1=(psi(i+1,j)+psi(i-1,j)+psi(i,j+1)+psi(i,j-1)+xi(i,j)*(dh^2))/4-psi(i,j);
            psi(i,j)=psi(i,j)+rho*err1;
            %vorticity
            err2=dt*((-dh/2)*(u(i,j)*(xi(i+1,j)-xi(i-1,j))+v(i,j)*(xi(i,j+1)-xi(i,j-1))) ...
                +(xi(i+1,j)+xi(i-1,j)+xi(i,j+1)+xi(i,j-1)-4*xi(i,j))/re)/(dh^2);
            xi(i,j)=xi(i,j)+rho*err2;
            temp=max(abs(err1),abs(err2));
            if err<temp
                err=temp;
            end
        end
    end

    if err<1e-4
        break;
    end
end


figure(1);
streamslice(v,u,10);
uv=zeros(n,n);
for i=1:n
    for j=1:n
        uv(i,j)=sqrt(u(i,j)^2+v(i,j)^2);
    end
end
figure(2);
contour(psi,100);
toc