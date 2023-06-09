clear all;clc;
tic
% gradE=[1,0,0;0,0,0;0,0,-1]; %[x; 0; -z]
dh=0.01;
[x,y,z] = meshgrid(-2:dh:2);
[xx,yy] = meshgrid(-2:dh:2);
n=length(x);
phi1 = 0.5*x.^2 + 0.5*y.^2 - z.^2;
phixy = 0.5*xx.^2 + 0.5*yy.^2;
phixz = 0.5*xx.^2 - yy.^2;



phi = phi1;
[ex,ey,ez] = gradient(-phi);
%quiver3(x,y,z,ex,ey,ez,3);axis equal;hold on;
axis equal;
contourf(xx,yy,phixy,30)




toc