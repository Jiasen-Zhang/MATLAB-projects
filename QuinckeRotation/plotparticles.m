% Plot
t0=200;
hv=[12,60];
for i=1:N
    xi=3*i-2;
    yi=3*i-1;
    zi=3*i;
    subplot(2,2,1);grid on;hold on;axis equal;
    plot3(datar(xi,:),datar(yi,:),datar(zi,:),'b');
    plot3(datar(1,(tmax-t0)/dt:tmax0),datar(2,(tmax-t0)/dt:tmax0),datar(3,(tmax-t0)/dt:tmax0),'r','LineWidth',1);
    scatter3(datar(xi,1),datar(yi,1),datar(zi,1),'*','r');
    scatter3(datar(xi,tmax0),datar(yi,tmax0),datar(zi,tmax0),'filled','b');
    xlabel('x');ylabel('y');zlabel('z');
    xlim([-2,2]);    ylim([-2,2]);    zlim([-2,2]);
    
    subplot(2,2,2);grid on;hold on;axis equal;
    %scatter(datar(xi,tmax0),datar(zi,tmax0),'filled','b');
    plot3(datar(xi,(tmax-t0)/dt:tmax0),datar(yi,(tmax-t0)/dt:tmax0),datar(zi,(tmax-t0)/dt:tmax0),'r');
    xlabel('x');ylabel('y');zlabel('z');
    xlim([-2,2]);    ylim([-2,2]);    zlim([-2,2]);

    subplot(2,2,3);grid on;hold on;axis equal;
    plot3(datar(xi,:),datar(yi,:),datar(zi,:),'b');
    plot3(datar(1,(tmax-t0)/dt:tmax0),datar(2,(tmax-t0)/dt:tmax0),datar(3,(tmax-t0)/dt:tmax0),'r','LineWidth',1);
    scatter3(datar(xi,1),datar(yi,1),datar(zi,1),'*','r');
    scatter3(datar(xi,tmax0),datar(yi,tmax0),datar(zi,tmax0),'filled','b');
    xlabel('x');ylabel('y');zlabel('z');
    view([0,0]);
    xlim([-2,2]);    ylim([-2,2]);    zlim([-2,2]);

    subplot(2,2,4);grid on;hold on;axis equal;
    %scatter(datar(xi,tmax0),datar(zi,tmax0),'filled','b');
    plot3(datar(xi,(tmax-t0)/dt:tmax0),datar(yi,(tmax-t0)/dt:tmax0),datar(zi,(tmax-t0)/dt:tmax0),'r');
    xlabel('x');ylabel('y');zlabel('z');
    view([0,0]);
    xlim([-2,2]);    ylim([-2,2]);    zlim([-2,2]);

end