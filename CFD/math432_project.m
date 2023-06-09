clear all;clc;
% math432 project by Jiasen Zhang

% x direction: v, j
% y direction: u, i
tic
% parameters
re=100;
Li=1;
Lj=1;
ni=50;%max(50,Li*round(sqrt(re)));
nj=50;%max(50,Lj*round(sqrt(re)));
hi=Li/ni;
hj=Lj/nj;
dt=min(2/re,1e-2);
err=1e-6;

% matrices
psi=zeros(ni,nj);
w=zeros(ni,nj);
Mi=zeros(ni,nj);
Mj=zeros(ni,nj);
u=zeros(ni,nj);
v=zeros(ni,nj);

% C-N
ri = dt/(2*re*hi*hi);
rj = dt/(2*re*hj*hj);
ei = ones(ni,1);
Ai0 = spdiags([ei -2*ei ei], [-1 0 1], ni, ni);
Ai1 = spdiags(ei,0,ni,ni) - ri * Ai0;
Ai2 = spdiags(ei,0,ni,ni) + ri * Ai0;
ej = ones(nj,1);
Aj0 = spdiags([ej -2*ej ej], [-1 0 1], nj, nj);
Aj1 = spdiags(ej,0,nj,nj) - rj * Aj0;
Aj2 = spdiags(ej,0,nj,nj) + rj * Aj0;

% Possion
IT = speye(nj-2);
eT = ones(ni-2,1);
T = spdiags([eT/(hi^2) -2*eT*((1/hi^2)+(1/hj^2)) eT/(hi^2)],[-1 0 1],ni-2,ni-2);
IS = speye(ni-2)/(hj^2);
eS = ones(nj-2,1);
S = spdiags([eS eS],[-1 1],nj-2,nj-2);
A = kron(IT,T) + kron(S,IS);


[x,y] = meshgrid(0:hj:Lj-hj,0:hi:Li-hi);
ii=2:ni-1;
ij=2:nj-1;
ii2=3:ni-2;
ij2=3:nj-2;
for k=1:1000000
    w0 = w; % old w
    
    % boundary
    w(:,1)=2*(psi(:,1)-psi(:,2))/(hj^2);
    w(:,nj)=2*(psi(:,nj)-psi(:,nj-1))/(hj^2);
    for j=1:nj
        %w(ni,j)=2*(psi(ni,j)-psi(ni-1,j))/(hi^2) + j*hi*(1-j*hi)*2/hi;
        w(ni,j)=2*(psi(ni,j)-psi(ni-1,j))/(hi^2) + 2/hi;
    end
    w(1,:)=2*(psi(1,:)-psi(2,:))/(hi^2);

    % corner
    w(1,1)=(w(1,2)+w(2,1))/2;
    w(1,nj)=(w(1,nj-1)+w(2,nj))/2;
    w(ni,1)=(w(ni,2)+w(ni-1,1))/2;
    w(ni,nj)=(w(ni,nj-1)+w(ni-1,nj))/2;
    
    %
    u(ii,ij)=(psi(ii,ij+1)-psi(ii,ij-1))/(2*hj);
    v(ii,ij)=-(psi(ii+1,ij)-psi(ii-1,ij))/(2*hi);
    uw = u.*w;
    vw = v.*w;
    
    %Lax-Wendroff
    Mi(ii,ij)=(uw(ii+1,ij)-uw(ii-1,ij))/(2*hi)...
        -dt*(uw(ii+1,ij)-2*uw(ii,ij)+uw(ii-1,ij))/(2*hi^2);
    Mj(ii,ij)=(vw(ii,ij+1)-vw(ii,ij-1))/(2*hj)...
        -dt*(vw(ii,ij+1)-2*vw(ii,ij)+vw(ii,ij-1))/(2*hj^2);
    
    % ADI
    w = Ai1\(w*Aj2-dt*(Mi+Mj)/2);   
    w = (Ai2*w-dt*(Mj+Mi)/2)/Aj1;

 	% LOD
%      w = Ai1\(Ai2*w-dt*Mi);  
%      w = (w*Aj2-dt*Mj)/Aj1; 
    
    % psi: interior points
    wvec = reshape(w(ii,ij),(ni-2)*(nj-2),1);
    psiv = -A\wvec;
    psi(ii,ij) = reshape(psiv,ni-2,nj-2);

    
    % stopping condition
    diff = max(max(abs(w-w0)));
    if mod(k,500)==0
        contour(x,y,psi);axis equal;drawnow;
        fprintf('%e, %e\n',err,diff);
    end
    if diff<err
        break;
    end
end

figure();
streamslice(x,y,v,u,10);axis equal;
figure();
contour(x,y,psi);axis equal;
figure();
subplot(1,2,1);
contour(x,y,u,'ShowText','on');axis equal;
title('v_y');
subplot(1,2,2);
contour(x,y,v,'ShowText','on');axis equal;
title('v_x');
toc