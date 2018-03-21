% Problem 16
clear all;
close all;
% set parametres
m=100;
n=20;
k=10000;
C=0.95;
tmax=40;
options=2;
nmax=1000;
alp=0.1;
eps=1e-8;
[~,~,~,~,~,~,Rall]=solver(m,n,k,C,tmax,options,nmax,alp,eps);
% plot Rall(R;t)
plot(Rall(2,:),Rall(1,:));
xlabel('time t')
ylabel('R(t)')
title('R(t) as a function of time from the initial condition up to steady state')