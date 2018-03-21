function [Y,u,v,phi,phi0,Rout,Rall]=solver(m,n,k,C,maxtime,options,nmax,alp,eps)
% set parameters
xe=5;
ye=1;
hx=xe/m;
hy=ye/n;
miu=2e-2;
alpha=1e-3;
% build x array for u xu and yu
xu=linspace(0,xe,m+1);
yu=linspace(-0.5*hy,ye+0.5*hy,n+2);
% build initial u array u0 no bc
u0nobc=initialu(xu,yu);
% add bc condition
u=ubc(u0nobc,yu);

% build initial v array v0 no bc
% build x array for u xu and yu
xv=linspace(-0.5*hx,xe+0.5*hx,m+2);
yv=linspace(0,ye,n+1);
% build initial v array v0 no bc
v0nobc=initialv(xv,yv);
% add bc condition
v=vbc(v0nobc,xv);

% build initial Y array
xy=linspace(0-2.5*hx,xe+2.5*hx,m+6);
yy=linspace(0-2.5*hy,ye+2.5*hy,n+6);
% build initial Y array v0 no bc
Y0nobc=initialY(xy,yy);
% apply bc to Y
Y=Ybc(Y0nobc,m,n,xy,yy);

% calculate R
R0=R(Y,m+6,n+6,hx,hy);
% initial time 
t=0;
% initial Huold and Hvold
Huold=0;
Hvold=0;
Rall(1,1)=R0;
% set initial guess phi
phi=zeros(m+2,n+2);
% time loop start
for q=1:k
% calculate dt
dtu=min(hx,hy)/(2*max(max(abs(u)))+max(max(abs(v))));
dtv=min(hx,hy)/(2*max(max(abs(v)))+max(max(abs(u))));
% dty= inifinity
dt=C*min(dtu,dtv);
Rall(2,q)=t;
% maxtime determination
if t<maxtime&&t+dt>maxtime
    dt=maxtime-t;
    Rall(2,q)=maxtime;
    t=t+dt
    % RK3/WENO5 + CN/ADI for Ynew
    [Y]=solveY2D(Y,u,v,alpha,hx,hy,dt,m,n,xy,yy);
    % Rnew
    Rout=R(Y,m+6,n+6,hx,hy);
    Rall(1,end)=Rout;
    break;
end

% steady state R determination
if q>20&&abs((Rall(1,q)-Rall(1,q-20))/(Rall(2,q)-Rall(2,q-20)))<eps
    % Rnew
    Rout=R(Y,m+6,n+6,hx,hy);
    Rall(1,end)=Rout;
    break;
end

% accumulate time t
t=t+dt

% RK3/WENO5 + CN/ADI for Ynew
[Y]=solveY2D(Y,u,v,alpha,hx,hy,dt,m,n,xy,yy);
% Rnew
Rnew=R(Y,m+6,n+6,hx,hy);
Rall(1,q+1)=Rnew;
% Adams-Bashforth/Crank-Nicholson for unew and vnew
[u,v,Huold,Hvold]=solveBurgers2D(u,v,hx,hy,dt,m,n,yu,xv,miu,q,Huold,Hvold);
% outlet correction of u (convective)
[u]=correction(u,v,hx,hy,m,n);
% rhs of possion equation as a source f(m+2,n+2)
[f]=possionrhs(u,v,hx,hy,dt,m,n);
% solve possion equation with f as a source phi(102x22)
[phi,~]=Vcycle(phi,f,m,hx,options,nmax,alp);
if q==1
    phi0=phi;
end
% project/correct velocities using Lagrange multiplier
[u,v]=projection(u,v,phi,dt,hx,hy,yu,xv,m,n);
end
end