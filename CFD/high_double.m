clear all;clc;
tic

%参
e=1e-6;
pr=1.5;
ra=100;
le=1;
lambda=1.4;
Li=1;Lj=2;
dh=1/20;  
ni=Li/dh;nj=Lj/dh;
dt=1e-4;

%变量设置          
u=zeros(ni,nj);
v=u;p=u;
u0=u;v0=v;
u1=u;v1=v;
u2=u;v2=v;
uj1=u;uj2=u;ui1=u;ui2=u;
vj1=v;vj2=v;vi1=v;vi2=v;
%变量设置 
s=zeros(ni,nj);t=s;
t0=t;s0=s;
t1=t;s1=s;
t2=t;s2=s;
tj1=t;tj2=t;ti1=t;ti2=t;
sj1=s;sj2=s;si1=s;si2=s;
D=zeros(ni,nj);        %散度

%初始条件
for i=1:ni
    for j=1:nj
        p(i,j)=1;
        s(i,j)=1-i*dh;
    end
end
%边界条件
for j=1:nj
    for i=1:ni
        s(i,1)=0;s(i,nj)=1;
        t(i,1)=0;t(i,nj)=1;
        %
        u(i,1)=0;
        v(1,j)=0;
        u(i,nj)=0;
        v(ni,j)=0;
        u(ni,j)=0;
        u(1,j)=0;
        v(i,1)=0;
        v(i,nj)=0;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:500000
    s0=s;t0=t;
    u0=u;v0=v;
    
    %求S T第一步
    sj1=compactd(s,1,1,dh);
    sj2=compactd(s,1,2,dh);
    si1=compactd(s,2,1,dh);
    si2=compactd(s,2,2,dh);
    tj1=compactd(t,1,1,dh);
    tj2=compactd(t,1,2,dh);
    ti1=compactd(t,2,1,dh);
    ti2=compactd(t,2,2,dh);
    for i=2:ni-1
        for j=2:nj-1
            errs=(sj2(i,j)+si2(i,j))/le-u(i,j)*si1(i,j)-v(i,j)*sj1(i,j);
            errt=tj2(i,j)+ti2(i,j)-u(i,j)*ti1(i,j)-v(i,j)*tj1(i,j);
            s(i,j)=errs*dt+s(i,j);
            t(i,j)=errt*dt+t(i,j);
        end
    end
  
    
    
    %求U V第一步
    vj1=compactd(v,1,1,dh);  %v j方向1阶导数
    vj2=compactd(v,1,2,dh);  %v j方向2阶导数
    vi1=compactd(v,2,1,dh);  %v i方向1阶导数
    vi2=compactd(v,2,2,dh);  %v i方向2阶导数
    uj1=compactd(u,1,1,dh);  %u j方向1阶导数
    uj2=compactd(u,1,2,dh);  %u j方向2阶导数
    ui1=compactd(u,2,1,dh);  %u i方向1阶导数
    ui2=compactd(u,2,2,dh);  %u i方向2阶导数
    for j=2:nj-1
        for i=2:ni-1
            fui=u(i,j)*ui1(i,j);
            fuj=v(i,j)*uj1(i,j);
            visu=(ui2(i,j)+uj2(i,j))*pr;
            fvi=u(i,j)*vi1(i,j);
            fvj=v(i,j)*vj1(i,j);
            visv=(vi2(i,j)+vj2(i,j))*pr;
            erru=((p(i,j)-p(i+1,j))/dh+visu-fui-fuj+ra*pr*(t(i,j)-lambda*s(i,j)));
            errv=((p(i,j)-p(i,j+1))/dh+visv-fvi-fvj);
            %
            u(i,j)=erru*dt+u(i,j);
            v(i,j)=errv*dt+v(i,j);
        end
    end
    
    %检验散度是否为零
    while(max(max((abs(D))))>e)        
        for j=2:nj
            for i=2:ni                         
                dp=-D(i,j)*dh*dh/(4*dt);
                dp2=dt*dp/dh;             
                u(i,j)=u(i,j)+dp2;
                u(i-1,j)=u(i-1,j)-dp2;
                v(i,j)=v(i,j)+dp2;
                v(i,j-1)=v(i,j-1)-dp2;
                D(i,j)=(u(i,j)-u(i-1,j)+v(i,j)-v(i,j-1))/(dh);                   
                p(i,j)=p(i,j)+dp;             
            end   
        end
    end
    for j=1:nj
        for i=1:ni
            s(i,1)=1;s(i,nj)=0;
            t(i,1)=1;t(i,nj)=0;
            %
            u(i,1)=0;
            v(1,j)=0;
            u(i,nj)=0;
            v(ni,j)=0;
            u(ni,j)=0;
            u(1,j)=0;
            v(i,1)=0;
            v(i,nj)=0;
        end
    end
    %检验收敛
    
    du=max(max(abs(u0-u)));dv=max(max(abs(v0-v)));
    ds=max(max(abs(s-s0)));dtt=max(max(abs(t-t0)));
    A=[du,dv,ds,dtt];
    err=max(A);
    if err<=e
        break;
    end
end
%计算完毕

%计算流函数
% w=zeros(n,n);
% psi=zeros(n,n);
% err1=w;
% for i=2:n-1
%     for j=2:n-1
%         w(i,j)=(v(i+1,j)-v(i-1,j))/(2*dh)-(u(i,j+1)-u(i,j-1))/(2*dh);
%     end
% end
% for k=1:100000
%     for i=2:n-1
%         for j=2:n-1
%             err1(i,j)=(psi(i+1,j)+psi(i-1,j)+psi(i,j+1)+psi(i,j-1)+w(i,j)*(dh^2))/4-psi(i,j);
%             psi(i,j)=psi(i,j)+err1(i,j);            
%         end
%     end
%     if max(max(abs(err1)))<1e-8
%         break;
%     end
% end
% contour(psi,n+1);hold on;

% uv=(u.^2+v.^2).^0.5;
% quiver(v./uv,u./uv,0.5);
streamslice(v,u,20);
axis equal;
axis tight;
toc