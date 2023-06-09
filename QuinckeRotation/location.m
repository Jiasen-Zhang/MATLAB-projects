fx=-15:1:15;
grid on;hold on;axis equal;
for i=1:N
     xi=3*i-2;
     yi=3*i-1;
     zi=3*i;
     scatter3(datar(xi,1),datar(yi,1),datar(zi,1),'o','r');
     scatter3(datar(xi,tmax0),datar(yi,tmax0),datar(zi,tmax0),'*','b');
     xlabel('x');ylabel('z');zlabel('y');
     legend('initial position','final position');
end
plot3(fx,0*fx,0*fx,'b');grid on;hold on;axis equal;
plot3(0*fx,fx,0*fx,'b');
plot3(0*fx,0*fx,fx,'b');view(-64,16);