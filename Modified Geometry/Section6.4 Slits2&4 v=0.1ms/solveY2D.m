function [Ynew]=solveY2D(Y0,u,v,alpha,hx,hy,dt,m,n,xy,yy)
% move u and v to the cell centers
u=reformu(u,m,n);
v=reformv(v,m,n);
% RK start
% Y0=Y0
% Y1
Y1=zeros(m+6,n+6);
% dYdx (no ghost cell)
dYdx0=dYdx(u,hx,Y0,m,n);
% dYdy
dYdy0=dYdy(v,hy,Y0,m,n);

for j=4:n+3
for i=4:m+3
Y1(i,j)=Y0(i,j)-1*(u(i-3,j-3)*dt*dYdx0(i-3,j-3))-1*(v(i-3,j-3)*dt*dYdy0(i-3,j-3));
end
end
% add ghost cell values
Y1=Ybc(Y1,m,n,xy,yy);

% Y2
Y2=zeros(m+6,n+6);
dYdx1=dYdx(u,hx,Y1,m,n);
dYdy1=dYdy(v,hy,Y1,m,n);

for j=4:n+3
for i=4:m+3
Y2(i,j)=Y1(i,j)+0.75*(u(i-3,j-3)*dt*dYdx0(i-3,j-3)+v(i-3,j-3)*dt*dYdy0(i-3,j-3))-0.25*(u(i-3,j-3)*dt*dYdx1(i-3,j-3)+v(i-3,j-3)*dt*dYdy1(i-3,j-3));
end
end
% add ghost cell values
Y2=Ybc(Y2,m,n,xy,yy);

% Y3
Y3=zeros(m+6,n+6);
dYdx2=dYdx(u,hx,Y2,m,n);
dYdy2=dYdy(v,hy,Y2,m,n);

for j=4:n+3
for i=4:m+3
Y3(i,j)=Y2(i,j)+1/12*(u(i-3,j-3)*dt*dYdx0(i-3,j-3)+v(i-3,j-3)*dt*dYdy0(i-3,j-3))+1/12*(u(i-3,j-3)*dt*dYdx1(i-3,j-3)+v(i-3,j-3)*dt*dYdy1(i-3,j-3))-2/3*(u(i-3,j-3)*dt*dYdx2(i-3,j-3)+v(i-3,j-3)*dt*dYdy2(i-3,j-3));
end
end
% add ghost cell values
Y3=Ybc(Y3,m,n,xy,yy);

%HY (HY 106x26) as a source to the parabolic solver
HY=(Y3-Y0)/dt;

% parabolic part
% get rid of bc of Ybc(100X20) and HY(100X20)
Y=Y0(3:end-2,3:end-2);
S=HY(3:end-2,3:end-2);
% set d1 and d2
d1=0.5*alpha*dt/hx^2;
d2=0.5*alpha*dt/hy^2;
% build A1,B1,C1 array when implicit x only
A1=-d1*ones(m,n);
B1=(1+2*d1)*ones(m,n);
C1=-d1*ones(m,n);
% build D1 array when implicit x only
D1=zeros(m,n);
for j=2:n+1
    for i=2:m+1
        D1(i-1,j-1)=d2*Y(i,j+1)+d2*Y(i,j-1)+(1-2*d2)*Y(i,j)+0.5*dt*S(i,j);
    end
end
% apply left boundary to A1,B1,C1,D1
A1(1,:)=0;
B1(1,:)=1+3*d1; 
C1(1,:)=-d1;
for j=2:n+1
    if yy(j+2)<0.5
D1(1,j-1)=d2*Y(2,j+1)+(1-2*d2)*Y(2,j)+d2*Y(2,j-1)+0.5*dt*S(2,j)+2*d1;
    else 
D1(1,j-1)=d2*Y(2,j+1)+(1-2*d2)*Y(2,j)+d2*Y(2,j-1)+0.5*dt*S(2,j);
    end
end
% apply right boundary to A1,B1,C1,D1
A1(end,:)=-d1;
B1(end,:)=1+d1;
C1(end,:)=0;
for j=2:n+1
D1(end,j-1)=d2*Y(end-1,j+1)+(1-2*d2)*Y(end-1,j)+d2*Y(end-1,j-1)+0.5*dt*S(end-1,j);
end
% calculate Y^(n+0.5) restore into Y
for j=1:n
    Y(2:end-1,j+1)=GE(A1(:,j),B1(:,j),C1(:,j),D1(:,j));
end
% add ghost cell values
Yhalf=zeros(m+6,n+6);
Yhalf(3:end-2,3:end-2)=Y;
Y=Ybc(Yhalf,m,n,xy,yy);

% step 2 initial
Y=Y(3:end-2,3:end-2);
% build A2,B2,C2 array when implicit x only
A2=-d2*ones(n,m);
B2=(1+2*d2)*ones(n,m);
C2=-d2*ones(n,m);
% build D2 array when implicit x only
D2=zeros(n,m);
for j=2:n+1
    for i=2:m+1
        D2(j-1,i-1)=d1*Y(i+1,j)+d1*Y(i-1,j)+(1-2*d1)*Y(i,j)+0.5*dt*S(i,j);
    end
end
% apply left boundary to A2,B2,C2,D2
A2(1,:)=0;
B2(1,:)=1+d2;
C2(1,:)=-d2;
for i=2:m+1
    if xy(i+2)>=0.875&&xy(i+2)<=1.125
            B2(1,i-1)=1+3*d2;
            D2(1,i-1)=d1*Y(i+1,2)+(1-2*d1)*Y(i,2)+d1*Y(i-1,2)+0.5*dt*S(i,2);
    else if xy(i+2)>=2.875&&xy(i+2)<=3.125
            B2(1,i-1)=1+3*d2;
            D2(1,i-1)=d1*Y(i+1,2)+(1-2*d1)*Y(i,2)+d1*Y(i-1,2)+0.5*dt*S(i,2);
             else if xy(i+2)>=1.875&&xy(i+2)<=2.125
            B2(1,i-1)=1+3*d2;
            D2(1,i-1)=d1*Y(i+1,2)+(1-2*d1)*Y(i,2)+d1*Y(i-1,2)+0.5*dt*S(i,2);
             else if xy(i+2)>=3.875&&xy(i+2)<=4.125
            B2(1,i-1)=1+3*d2;
            D2(1,i-1)=d1*Y(i+1,2)+(1-2*d1)*Y(i,2)+d1*Y(i-1,2)+0.5*dt*S(i,2);
                 end
                 end
        end
    end
end
% apply right boundary to A2,B2,C2,D2
A2(end,:)=-d2;
B2(end,:)=1+d2;
C2(end,:)=0;
for i=2:m+1
    if xy(i+2)>=1.875&&xy(i+2)<=2.125
                B2(end,i-1)=1+3*d2;
                D2(end,i-1)=d1*Y(i+1,end-1)+(1-2*d1)*Y(i,end-1)+d1*Y(i-1,end-1)+0.5*dt*S(i,end-1)+2*d2;
    else if xy(i+2)>=3.875&&xy(i+2)<=4.125
                B2(end,i-1)=1+3*d2;
                D2(end,i-1)=d1*Y(i+1,end-1)+(1-2*d1)*Y(i,end-1)+d1*Y(i-1,end-1)+0.5*dt*S(i,end-1)+2*d2;
                 else if xy(i+2)>=0.875&&xy(i+2)<=1.125
                B2(end,i-1)=1+3*d2;
                D2(end,i-1)=d1*Y(i+1,end-1)+(1-2*d1)*Y(i,end-1)+d1*Y(i-1,end-1)+0.5*dt*S(i,end-1)+2*d2;
                 else if xy(i+2)>=2.875&&xy(i+2)<=3.125
                B2(end,i-1)=1+3*d2;
                D2(end,i-1)=d1*Y(i+1,end-1)+(1-2*d1)*Y(i,end-1)+d1*Y(i-1,end-1)+0.5*dt*S(i,end-1)+2*d2;
                     end
                     end
        end
    end
end
% calculate Y^(n+1) restore into Y
for i=1:m
    Y(i+1,2:end-1)=GE(A2(:,i),B2(:,i),C2(:,i),D2(:,i));
end
% add bc values  to Y
Ynew=zeros(m+6,n+6);
Ynew(3:end-2,3:end-2)=Y;
Ynew=Ybc(Ynew,m,n,xy,yy);
end