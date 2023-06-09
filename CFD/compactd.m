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

