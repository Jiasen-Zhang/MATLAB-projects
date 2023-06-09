clear all;clc;
tic
% gradE=[1,0,0;0,0,0;0,0,-1]; %[x; 0; -z]
dh=0.1;
[x,z] = meshgrid(-2:dh:2);
n=length(x);
phi1 = x.^2 - z.^2;
phi2 = x.^3 - 3*x.*z.^2;
phi3 = x.^4 - 6*x.^2.*z.^2+z.^4;

phi = phi1;
[ex,ez] = gradient(-phi);
%streamslice(x,z,ex,ez);hold on;
contourf(x,z,phi,50);





toc