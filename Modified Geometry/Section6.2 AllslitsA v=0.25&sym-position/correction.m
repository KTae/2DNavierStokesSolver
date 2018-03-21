% Outlet Correction
function [u]=correction(u,v,dx,dy,m,n)
ly=1;
A=0;
B=0;
C=0;
D=0;
for j=1:n
A=A+(max(0,u(1,j+1))+max(0,-u(m+1,j+1)))*dy;
C=C+(max(0,-u(1,j+1))+max(0,u(m+1,j+1)))*dy;
end
for i=1:m
B=B+(max(0,v(i+1,1))+max(0,-v(i+1,n+1)))*dx;
D=D+(max(0,-v(i+1,1))+max(0,v(i+1,n+1)))*dx;
end
qin=A+B;
qout=C+D;
qcorr=qin-qout;
ucorr=qcorr/ly;
u(m+1,:)=u(m+1,:)+ucorr;
end