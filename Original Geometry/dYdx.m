% calculate dY/dx
function w=dYdx(u,hx,Y,m,n)

A=zeros(m,n);
B=zeros(m,n);
C=zeros(m,n);
D=zeros(m,n);
w=zeros(m,n);

for j=4:n+3
for i=4:m+3
% u>0
if u(i-3,j-3)>0
    A(i-3,j-3)=1/hx*((Y(i-1,j)-Y(i-2,j))-(Y(i-2,j)-Y(i-3,j)));
    B(i-3,j-3)=1/hx*((Y(i,j)-Y(i-1,j))-(Y(i-1,j)-Y(i-2,j)));
    C(i-3,j-3)=1/hx*((Y(i+1,j)-Y(i,j))-(Y(i,j)-Y(i-1,j)));
    D(i-3,j-3)=1/hx*((Y(i+2,j)-Y(i+1,j))-(Y(i+1,j)-Y(i,j)));
    w(i-3,j-3)=1/(12*hx)*(-(Y(i-1,j)-Y(i-2,j))+7*(Y(i,j)-Y(i-1,j))+7*(Y(i+1,j)-Y(i,j))-(Y(i+2,j)-Y(i+1,j)))-weno(A(i-3,j-3),B(i-3,j-3),C(i-3,j-3),D(i-3,j-3));
end
% u<0
if u(i-3,j-3)<0
    A(i-3,j-3)=1/hx*((Y(i+3,j)-Y(i+2,j))-(Y(i+2,j)-Y(i+1,j)));
    B(i-3,j-3)=1/hx*((Y(i+2,j)-Y(i+1,j))-(Y(i+1,j)-Y(i,j)));
    C(i-3,j-3)=1/hx*((Y(i+1,j)-Y(i,j))-(Y(i,j)-Y(i-1,j)));
    D(i-3,j-3)=1/hx*((Y(i,j)-Y(i-1,j))-(Y(i-1,j)-Y(i-2,j)));
    w(i-3,j-3)=1/(12*hx)*(-(Y(i-1,j)-Y(i-2,j))+7*(Y(i,j)-Y(i-1,j))+7*(Y(i+1,j)-Y(i,j))-(Y(i+2,j)-Y(i+1,j)))+weno(A(i-3,j-3),B(i-3,j-3),C(i-3,j-3),D(i-3,j-3));
end
end
end
end
