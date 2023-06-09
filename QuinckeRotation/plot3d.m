t0=100;
for i=1:N
     xi=3*i-2;
     yi=3*i-1;
     zi=3*i;
     subplot(2,2,1);grid on;hold on;axis equal;
     plot3(datar(xi,:),datar(yi,:),datar(zi,:));grid on;hold on;axis equal;
     %plot3(datar(1,(tmax-t0)/dt:tmax0),datar(2,(tmax-t0)/dt:tmax0),datar(3,(tmax-t0)/dt:tmax0),'LineWidth',1);
     scatter3(datar(xi,1),datar(yi,1),datar(zi,1),'*');
     scatter3(datar(xi,tmax0),datar(yi,tmax0),datar(zi,tmax0),'filled');
     xlabel('x');ylabel('y');zlabel('z');
     % 
%      plot(datar(xi,:),datar(zi,:));grid on;hold on;axis equal;
%      plot(datar(1,(tmax-t0)/dt:tmax0),datar(3,(tmax-t0)/dt:tmax0),'LineWidth',1);
%      scatter(datar(xi,1),datar(zi,1),'*');
%      scatter(datar(xi,tmax0),datar(zi,tmax0),'filled');
%      xlabel('x');ylabel('z');

     subplot(2,2,2);grid on;hold on;
     plot(time,dataom3(i,:));
     xlabel('t');ylabel('\Omega');
     
     subplot(2,2,3);grid on;hold on;
     plot(time,datar0(i,:));
     xlabel('t');ylabel('distance');
     
     subplot(2,2,4);axis equal;grid on;hold on;
     scatter3(datar(xi,tmax0),datar(yi,tmax0),datar(zi,tmax0),'filled');
     plot3(datar(xi,(tmax-t0)/dt:tmax0),datar(yi,(tmax-t0)/dt:tmax0),datar(zi,(tmax-t0)/dt:tmax0),'LineWidth',1);
     xlabel('x');ylabel('y');zlabel('z');

end