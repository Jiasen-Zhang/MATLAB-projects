function yp=F1(t,y)
yp=zeros(12,1);
r=y(1:6,1);
P=y(7:12,1);
Om=fom1(r,P);
yp(1:6,1)=fr1(r,P,Om);
yp(7:12,1)=fp1(r,P,Om);

end
