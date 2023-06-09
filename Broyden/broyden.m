function [ bestsolution,fmin,k ] =broyden(f,phi,x0,e)
%broyden族函数,指定参数phi的值即可
%f 目标函数
%x0 初始点
%e 允许误差
%lamda 步长
%phi 加权组合参数
%bestsolution 最优解
%fmin 函数极小值
%k 迭代次数

h0=[1,0;0,1];
g0=mygradient(f,x0);%梯度函数
k=1;
while(1)
    dk=-h0*g0;
    f0=@(lamda)f(x0+lamda*dk);
    lamdak=goldensection(f0);%黄金分割法函数
    xk=x0+lamdak*dk;
    if abs(norm(mygradient(f,xk)))<=e
        bestsolution=xk;
        fmin=f(xk);     
        break;
    else
        if k==10
            x0=xk;
            h0=[1,0;0,1];
            g0=mygradient(f,x0);
            k=1;
        else
        gk=mygradient(f,xk);
        pk=xk-x0;
        qk=mygradient(f,xk)-mygradient(f,x0);
        vk=sqrt(qk'*h0*qk)*((pk/(pk'*qk))-((h0*qk)/(qk'*h0*qk)));
        hk=h0+((pk*pk')/(pk'*qk))-((h0*qk*qk'*h0)/(qk'*h0*qk))+(phi*vk*vk');
        x0=xk;
        g0=gk;
        h0=hk;   
        k=k+1;
        end
     
    end
    
end

end



