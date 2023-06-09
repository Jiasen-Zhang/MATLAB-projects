function [dp,dq] = fpq1(P,Q,r,Om)
%
global N;
global ecm;
global ecm1;
global scm;
global scm1;
global D;
global D1;
global G;
%global gradE;
%
dp=zeros(3*N,1);
dq=zeros(3*N,3);
Ep=zeros(3*N,N);
Eq=zeros(3*N,N);
Gp=zeros(3*N,3,N);
% Get Epij, Eqij, Gpij
for i=1:N
    for j=[1:i-1,i+1:N]% without j=i
        xi=3*i-2;
        xj=3*j-2;
        zi=3*i;
        zj=3*j;
        R=r(xj:zj)-r(xi:zi); %3x1 Rij=rj-ri
        normR=norm(R);%1x1
        R0=R/normR;%3x1
        Pj=P(xj:zj);%3x1
        Qj=Q(xj:zj,:);%3x3
        %
        Ep(xi:zi,j)=  (Pj - 3*R0'*Pj*R0)/(normR^3)  ;
        Eq(xi:zi,j)=  (Qj*R0 - 2.5*(R0'*Qj*R0)*R0)/(normR^4);
        Gpij=zeros(3,3);
        for ii=1:3
            Gpij(ii,ii)=(   -6*R(ii)*Pj(ii)/(normR^5) - 3*R'*Pj/(normR^5) + (15*R'*Pj*R(ii)^2)/(normR^7)   );
            for jj=[1:ii-1,ii+1:3]
                Gpij(ii,jj)=(   -3*(R(ii)*Pj(jj)+R(jj)*Pj(ii))/(normR^5) + (15*R'*Pj*R(ii)*R(jj))/(normR^7)   );
            end
        end
        Gp(xi:zi,:,j)=Gpij;
    end
end
sumE=sum(Ep+Eq,2);
sumGp=0*sum(Gp,3);
% Get dp and dq
for i=1:N
    xi=3*i-2;
    yi=3*i-1;
    zi=3*i;
    gradE=[0,0,1;
        0,0,0;
        1,0,0];
    Ei=[r(zi);0;r(xi)];

    %
    dPi1=P(xi:zi)+ecm*(-G*Ei+sumE(xi:zi));
    dPi2=(P(xi:zi)+scm*(-G*Ei+sumE(xi:zi)))/D;
    dp(xi:zi)=cross(Om(xi:zi),dPi1)-dPi2;
    %
    dQi1=Q(xi:zi,:)+2*ecm1*(-G*gradE+sumGp(xi:zi,:));
    dQi2=[cross(Om(xi:zi),dQi1(1,:));cross(Om(xi:zi),dQi1(2,:));cross(Om(xi:zi),dQi1(3,:))];
    dQi3=(Q(xi:zi,:)+2*scm1*(-G*gradE+sumGp(xi:zi,:)))/D1;
    dq(xi:zi,:)=dQi2+dQi2'-dQi3;
end
end