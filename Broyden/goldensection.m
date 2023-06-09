
function lamdak = goldensection(f)
%���ý��˷�Ѱ�ҳ�ʼ����
%x1 ������ʼ��
%h ��ʼ����
x1=0;
h=1;
f1=f(x1);
k=0;
while 1
    x4=x1+h;
    f4=f(x4);
    k=k+1;
    if f4<f1
        x2=x1;x1=x4;f2=f1;f1=f4;
        h=2*h;
    else
        if k==1
            h=-h;x2=x4;f2=f4;
        else
            x3=x2;x2=x1;x1=x4;break;
        end
    end
end
%a,b ��ʼ����
if x3>x1
    a=x1;b=x3;
else
    a=x3;b=x1;
end
%����������lamdak�Ļƽ�ָ
%l �������
l=0.001;
for k=1:100
    m=a+0.382*(b-a);
    n=a+0.618*(b-a);
    if f(m)>f(n)
        a=m;m=n;
    else b=n;n=m;
    end
    if abs(b-a)<l
       lamdak=0.5*(a+b);
        break;
    end
end

