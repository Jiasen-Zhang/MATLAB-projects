function main
%13307110108 张加森
%该段程序用于计算不同phi下的最优解，相应函数值和迭代次数
%并画出求得函数极小值与参数phi的散点图
clc;
i=1;
f=@(x)0.5*x(1)^2+x(2)^2-x(1)*x(2);%设置函数表达式
x0=[2;2];
e=1e-05;
%定义三个矩阵
best=eye(2,10);
fmini=eye(1,10);
ki=eye(1,10);
for phi=0:0.1:1 %phi从0到1，间隔为0.1
    [ bestsolution,fmin,k ] =broyden(f,phi,x0,e);
    best(:,i)=bestsolution;
    fmini(i)=fmin;
    ki(i)=k;
    scatter(phi,fmin,10,'filled');%画出极小值fmin与phi的关系散点图
    hold on;
    i=i+1;
end
best%依次列出各phi值下的最优解
fmini%依次列出各phi值下的fmin
ki%依次列出各phi值下的迭代次数
end

