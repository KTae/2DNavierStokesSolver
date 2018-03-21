% Problem 15 
% set parameters
clear all;
close all;
m=100;
n=20;
k=100000;
C=0.95;
eps=1e-7;
tmax=[0.05 0.1 1 2 10] ;
options=2;
nmax=1500;
alp=0.1;
xe=5;
ye=1;
hx=xe/m;
hy=ye/n;
% run the solver
Y=zeros(m+6,n+6,length(tmax)+1);
u=zeros(m+1,n+2,length(tmax)+1);
v=zeros(m+2,n+1,length(tmax)+1);
phi=zeros(m+2,n+2,length(tmax)+1);
for i=1:length(tmax)
[Y(:,:,i+1),u(:,:,i+1),v(:,:,i+1),phi(:,:,i+1),~,~]=solver(m,n,k,C,tmax(i),options,nmax,alp,eps);
end
% restore initial guess Y,U,V,PHI
% build x array for u xu and yu
xu=linspace(0,xe,m+1);
yu=linspace(-0.5*hy,ye+0.5*hy,n+2);
% build initial u array u0 no bc
u0nobc=initialu(xu,yu);
% add bc condition
u(:,:,1)=ubc(u0nobc,yu);

% build initial v array v0 no bc
% build x array for u xu and yu
xv=linspace(-0.5*hx,xe+0.5*hx,m+2);
yv=linspace(0,ye,n+1);
% build initial v array v0 no bc
v0nobc=initialv(xv,yv);
% add bc condition
v(:,:,1)=vbc(v0nobc,xv);

% build initial Y array
xy=linspace(0-2.5*hx,xe+2.5*hx,m+6);
yy=linspace(0-2.5*hy,ye+2.5*hy,n+6);
% build initial Y array v0 no bc
Y0nobc=initialY(xy,yy);
% apply bc to Y
Y(:,:,1)=Ybc(Y0nobc,m,n,xy,yy);

% build initial phi array(m+2,n+2)
phix=xy(3:end-2);
phiy=yy(3:end-2);
phi(:,:,1)=zeros(m+2,n+2);

% Y(1-Y) in Q
Q=zeros(m+6,n+6,5);
for i=1:length(tmax)+1
Q(:,:,i)=Y(:,:,i).*(1-Y(:,:,i));
end

% plot figure
% plot u
figure(1)
for i=1:length(tmax)+1
subplot(3,2,i)
contourf(xu,yu,u(:,:,i)');
axis([0 5 0 1]);
colorbar
xlabel('x coordinate[cm]')
ylabel('y coordinate[cm]')
title('U[m/s] time=  s')
end
% plot v
figure(2)
for i=1:length(tmax)+1
subplot(3,2,i)
contourf(xv,yv,v(:,:,i)');
axis([0 5 0 1]);
colorbar
xlabel('x coordinate[cm]')
ylabel('y coordinate[cm]')
title('V[m/s] time=  s')
end
% plot Y
figure(3)
for i=1:length(tmax)+1
subplot(3,2,i)
contourf(xy,yy,Y(:,:,i)');
axis([0 5 0 1]);
colorbar
xlabel('x coordinate[cm]')
ylabel('y coordinate[cm]')
title('Y[-] time=  s')
end
% Plot Q
figure(4)
for i=1:length(tmax)+1
subplot(3,2,i)
contourf(xy,yy,Q(:,:,i)');
axis([0 5 0 1]);
colorbar
xlabel('x coordinate[cm]')
ylabel('y coordinate[cm]')
title('Y(1-Y)[-] time=  s')
end
% plot phi
figure(5)
for i=1:length(tmax)+1
subplot(3,2,i)
contourf(phix,phiy,phi(:,:,i)');
axis([0 5 0 1]);
colorbar
xlabel('x coordinate[cm]')
ylabel('y coordinate[cm]')
title('Phi[-] time=  s')
end