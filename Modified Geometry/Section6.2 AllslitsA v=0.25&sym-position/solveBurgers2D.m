function [unew,vnew,Huold,Hvold]=solveBurgers2D(u,v,dx,dy,dt,m,n,yu,xv,miu,q,Huold,Hvold)
% for u
duudx=zeros(m+1,n+2);
duvdy=zeros(m+1,n+2);
for j=2:n+1
    for i=2:m
        duudx(i,j)=((u(i+1,j)+u(i,j))^2-(u(i,j)+u(i-1,j))^2)/(4*dx);
        duvdy(i,j)=((u(i,j)+u(i,j+1))*(v(i+1,j)+v(i,j))-(u(i,j-1)+u(i,j))*(v(i,j-1)+v(i+1,j-1)))/(4*dy);
    end
end
% add boundary
duudx=ubc(duudx,yu);
duvdy=ubc(duvdy,yu);
% Hu
Hu=zeros(m+1,n+2);
for j=1:n+2
    for i=1:m+1
     Hu(i,j)=-duudx(i,j)-duvdy(i,j);
    end
end

% for v
dvvdy=zeros(m+2,n+1);
duvdx=zeros(m+2,n+1);
for j=2:n
    for i=2:m+1
        dvvdy(i,j)=((v(i,j+1)+v(i,j))^2-(v(i,j)+v(i,j-1))^2)/(4*dy);
        duvdx(i,j)=((u(i,j)+u(i,j+1))*(v(i+1,j)+v(i,j))-(u(i-1,j)+u(i-1,j+1))*(v(i-1,j)+v(i,j)))/(4*dx);
    end
end
% add boundary
dvvdy=vbc(dvvdy,xv);
duvdx=vbc(duvdx,xv);
% Hv
Hv=zeros(m+2,n+1);
for j=1:n+1
    for i=1:m+2
     Hv(i,j)=-dvvdy(i,j)-duvdx(i,j);
    end
end

% Huold and Hvold
if q==1
Huold=Hu;
Hvold=Hv;
end
% Source array Su and Sv
Su=1.5*Hu-0.5*Huold;
Sv=1.5*Hv-0.5*Hvold;

% hyperbolic part
% set d1 and d2
d1=0.5*miu*dt/dx^2;
d2=0.5*miu*dt/dy^2;

% ADI for u
% u step1
% build Au1,Bu1,Cu1 array when implicit x only
Au1=-d1*ones(m-1,n);
Bu1=(1+2*d1)*ones(m-1,n);
Cu1=-d1*ones(m-1,n);
% build Du1 array when implicit x only
Du1=zeros(m-1,n);
for j=2:n+1
    for i=2:m
        Du1(i-1,j-1)=d2*u(i,j+1)+d2*u(i,j-1)+(1-2*d2)*u(i,j)+0.5*dt*Su(i,j);
    end
end    
% apply left boundary to Au1,Bu1,Cu1,Du1
Au1(1,:)=0;
Bu1(1,:)=1+2*d1; 
Cu1(1,:)=-d1;
for j=2:n+1
Du1(1,j-1)=d2*u(2,j+1)+(1-2*d2)*u(2,j)+d2*u(2,j-1)+0.5*dt*Su(2,j)+d1*u(1,j);
end
% apply right boundary to Au1,Bu1,Cu1,Du1
Au1(end,:)=-d1;
Bu1(end,:)=1+d1;
Cu1(end,:)=0;
for j=2:n+1
Du1(end,j-1)=d2*u(end-1,j+1)+(1-2*d2)*u(end-1,j)+d2*u(end-1,j-1)+0.5*dt*Su(end-1,j);
end
% calculate u^(n+0.5) restore into u
for j=1:n
    u(2:end-1,j+1)=GE(Au1(:,j),Bu1(:,j),Cu1(:,j),Du1(:,j));
end
% add bc to u^(n+0.5)
us=ubc(u,yu);

% u step 2
% build Au2,Bu2,Cu2 array when implicit x only
Au2=-d2*ones(n,m-1);
Bu2=(1+2*d2)*ones(n,m-1);
Cu2=-d2*ones(n,m-1);
% build Du2 array when implicit x only
Du2=zeros(n,m-1);
for j=2:n+1
    for i=2:m
        Du2(j-1,i-1)=d1*us(i+1,j)+d1*us(i-1,j)+(1-2*d1)*us(i,j)+0.5*dt*Su(i,j);
    end
end
% apply left boundary to Au2,Bu2,Cu2,Du2
Au2(1,:)=0;
Bu2(1,:)=1+3*d2;
Cu2(1,:)=-d2;

% apply right boundary to Au2,Bu2,Cu2,Du2
Au2(end,:)=-d2;
Bu2(end,:)=1+3*d2;
Cu2(end,:)=0;

% calculate u^(n+1) restore into u
unew=zeros(m+1,n+2);
for i=1:m-1
    unew(i+1,2:end-1)=GE(Au2(:,i),Bu2(:,i),Cu2(:,i),Du2(:,i));
end
% add bc to unew
unew=ubc(unew,yu);

% ADI for v(m+2,n+1)
% v step1
% build Av1,Bv1,Cv1 array when implicit x only
Av1=-d1*ones(m,n-1);
Bv1=(1+2*d1)*ones(m,n-1);
Cv1=-d1*ones(m,n-1);
% build Dv1 array when implicit x only
Dv1=zeros(m,n-1);
for j=2:n
    for i=2:m+1
        Dv1(i-1,j-1)=d2*v(i,j+1)+d2*v(i,j-1)+(1-2*d2)*v(i,j)+0.5*dt*Sv(i,j);
    end
end    
% apply left boundary to Av1,Bv1,Cv1,Dv1
Av1(1,:)=0;
Bv1(1,:)=1+3*d1; 
Cv1(1,:)=-d1;

% apply right boundary to Av1,Bv1,Cv1,Dv1
Av1(end,:)=-d1;
Bv1(end,:)=1+3*d1;
Cv1(end,:)=0;

% calculate v^(n+0.5) restore into vs
for j=1:n-1
    v(2:end-1,j+1)=GE(Av1(:,j),Bv1(:,j),Cv1(:,j),Dv1(:,j));
end
% add bc to v^(n+0.5) vs
vs=vbc(v,xv);

% v step2
% build Av2,Bv2,Cv2 array when implicit x only
Av2=-d2*ones(n-1,m);
Bv2=(1+2*d2)*ones(n-1,m);
Cv2=-d2*ones(n-1,m);
% build Dv2 array when implicit x only
Dv2=zeros(n-1,m);
for j=2:n
    for i=2:m+1
        Dv2(j-1,i-1)=d1*vs(i+1,j)+d1*vs(i-1,j)+(1-2*d1)*vs(i,j)+0.5*dt*Sv(i,j);
    end
end

% apply left boundary to Av2,Bv2,Cv2,Dv2
Av2(1,:)=0;
Bv2(1,:)=1+2*d2;
Cv2(1,:)=-d2;
for i=2:m+1
    if xv(i)>=0.875&&xv(i)<=1.125
            Dv2(1,i-1)=d1*vs(i+1,2)+(1-2*d1)*vs(i,2)+d1*vs(i-1,2)+0.5*dt*Sv(i,2)+0.25*d2;
    else if xv(i)>=2.875&&xv(i)<=3.125
            Dv2(1,i-1)=d1*vs(i+1,2)+(1-2*d1)*vs(i,2)+d1*vs(i-1,2)+0.5*dt*Sv(i,2)+0.25*d2;
        end
    end
end
% apply right boundary to A2,B2,C2,D2
Av2(end,:)=-d2;
Bv2(end,:)=1+2*d2;
Cv2(end,:)=0;
for i=2:m+1
    if xv(i)>=0.875&&xv(i)<=1.125
                Dv2(end,i-1)=d1*vs(i+1,end-1)+(1-2*d1)*vs(i,end-1)+d1*vs(i-1,end-1)+0.5*dt*Sv(i,end-1)-0.25*d2;
    else if xv(i)>=2.875&&xv(i)<=3.125
                Dv2(end,i-1)=d1*vs(i+1,end-1)+(1-2*d1)*vs(i,end-1)+d1*vs(i-1,end-1)+0.5*dt*Sv(i,end-1)-0.25*d2;
        end
    end
end

% calculate v^(n+1) restore into vs
for j=1:m
    vs(j+1,2:end-1)=GE(Av2(:,j),Bv2(:,j),Cv2(:,j),Dv2(:,j));
end
% add bc to v^(n+1) restore in vnew
vnew=vbc(vs,xv);
end