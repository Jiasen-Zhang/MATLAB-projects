function Om = fomn(P,Q,r)
global N;
global G;
%
Om=zeros(3*N,1);
Ep=zeros(3*N,N);
Eq=zeros(3*N,N);
Gp=zeros(3*N,3,N);
Fd2i0=zeros(3*N,N);
dOmi1=zeros(3*N,N);
Fel=zeros(3*N,1);
Tel=zeros(3*N,1);
%
for i=1:N
    for j=[1:i-1,i+1:N]% without j=i
        xi=3*i-2;
        xj=3*j-2;
        zi=3*i;
        zj=3*j;
        R=r(xj:zj)-r(xi:zi); %3x1 Rij=rj-ri
        normR=norm(R);%1x1
        R0=R/normR;%3x1
        Pj=P(xj:zj); Pi=P(xi:zi);%3x1
        Qj=Q(xj:zj,:);%3x3
        %
        Fd2i0m= Pi'*R0*Pj + Pj'*R0*Pi + Pi'*Pj*R0 - 5*(Pj'*R0)*(Pi'*R0)*R0;
        Fd2i0(xi:zi,j)=((-1)^i)*12*pi*(Fd2i0m)/(normR^4);
        %
        Ep(xi:zi,j)=(  Pj/(normR^3) - 3*R'*Pj*R/(normR^5)  );
        Eq(xi:zi,j)=(  Qj*R/(normR^5) - 2.5*(R'*Qj*R)*R/(normR^7)  );
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
sumFd2i0=sum(Fd2i0,2);
sumEp=sum(Ep,2);
sumEq=sum(Eq,2);
sumGp=0*sum(Gp,3);
% Get Fel and Tel
for i=1:N
    xi=3*i-2;
    yi=3*i-1;
    zi=3*i;
    
    gradE=[1,0,0;
        0,0,0;
        0,0,-1];
    Ei=[r(xi);
        0;
        -r(zi)];

    %
    Fd1i=-4*pi*G*(-gradE)*P(xi:zi);
    Fd2i=-sumFd2i0(xi:zi);
    Fd3i=0;
    Fel(xi:zi)=Fd1i+Fd2i+Fd3i;
    %
    Td1i=-4*pi*G*(cross(P(xi:zi),-Ei)-[Q(yi,:)*gradE(3,:)'-Q(zi,:)*gradE(2,:)';
        Q(zi,:)*gradE(1,:)'-Q(xi,:)*gradE(3,:)';
        Q(xi,:)*gradE(2,:)'-Q(yi,:)*gradE(1,:)']);
    Td2i=-4*pi*cross(P(xi:zi),sumEp(xi:zi));
    Td3i=-4*pi*cross(P(xi:zi),sumEq(xi:zi));
    Td4i=-4*pi*[Q(yi,:)*sumGp(zi,:)'-Q(zi,:)*sumGp(yi,:)';
        Q(zi,:)*sumGp(xi,:)'-Q(xi,:)*sumGp(zi,:)';
        Q(xi,:)*sumGp(yi,:)'-Q(yi,:)*sumGp(xi,:)'];
    Tel(xi:zi)=Td1i+Td2i+Td3i+Td4i;
end
% Get dui1, dui2, dOmi1
for i=1:N
    for j=[1:i-1,i+1:N]% without j=i
        xi=3*i-2;
        xj=3*j-2;
        zi=3*i;
        zj=3*j;
        R=r(xj:zj)-r(xi:zi); %3x1 Rij=rj-ri
        normR=norm(R);%1x1
        R0=R/normR;%3x1
        %
        dOmi1(xj:zj,j)= -Tel(xj:zj)/(16*pi*normR^3) - 3*Tel(xj:zj)'*R0*R0/(16*pi*normR^3) - cross(Fel(xj:zj),R0)/(8*pi*normR^2);
    end
end
sumdOmi1=sum(dOmi1,2);
% Get ri, Omi
for i=1:N
    xi=3*i-2;
    zi=3*i;
    Om(xi:zi)=Tel(xi:zi)/(8*pi)+sumdOmi1(xi:zi);
end

end