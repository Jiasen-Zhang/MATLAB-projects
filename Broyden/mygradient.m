
function g = mygradient( f,x0 )
%��ַ��ݶȺ���
%ȡ����1e-05
delta=1e-05;
%�ֱ���
g1=(f(x0+[delta;0])-f(x0-[delta;0]))/(2*delta);
g2=(f(x0+[0;delta])-f(x0-[0;delta]))/(2*delta);
g=[g1;g2];
end

