function output=F(t,data)
global N;
r = data(1:3*N,1);
P = data(3*N+1:6*N,1);
Q(:,1) = data(6*N+1:9*N,1);
Q(:,2) = data(9*N+1:12*N,1);
Q(:,3) = data(12*N+1:15*N,1);

output = zeros(18*N,1);

Om = fomn(P,Q,r);
dr=frn(P,Q,r);
[dp,dq]=fpqn(P,Q,r,Om);

output(1:3*N,1)=dr;
output(3*N+1:6*N,1)=dp;
output(6*N+1:9*N)=dq(:,1);
output(9*N+1:12*N)=dq(:,2);
output(12*N+1:15*N)=dq(:,3);
output(15*N+1:18*N,1) = Om;

end