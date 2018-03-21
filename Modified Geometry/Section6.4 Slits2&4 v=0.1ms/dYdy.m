% calculate dY/dy
function w=dYdy(v,hy,Y,m,n)
A=zeros(m,n);
B=zeros(m,n);
C=zeros(m,n);
D=zeros(m,n);
w=zeros(m,n);

for j=4:n+3
for i=4:m+3
% v>0
if v(i-3,j-3)>0
    A(i-3,j-3)=1/hy*((Y(i,j-1)-Y(i,j-2))-(Y(i,j-2)-Y(i,j-3)));
    B(i-3,j-3)=1/hy*((Y(i,j)-Y(i,j-1))-(Y(i,j-1)-Y(i,j-2)));
    C(i-3,j-3)=1/hy*((Y(i,j+1)-Y(i,j))-(Y(i,j)-Y(i,j-1)));
    D(i-3,j-3)=1/hy*((Y(i,j+2)-Y(i,j+1))-(Y(i,j+1)-Y(i,j)));
    w(i-3,j-3)=1/(12*hy)*(-(Y(i,j-1)-Y(i,j-2))+7*(Y(i,j)-Y(i,j-1))+7*(Y(i,j+1)-Y(i,j))-(Y(i,j+2)-Y(i,j+1)))-weno(A(i-3,j-3),B(i-3,j-3),C(i-3,j-3),D(i-3,j-3));
end
% v<0
if v(i-3,j-3)<0
    A(i-3,j-3)=1/hy*((Y(i,j+3)-Y(i,j+2))-(Y(i,j+2)-Y(i,j+1)));
    B(i-3,j-3)=1/hy*((Y(i,j+2)-Y(i,j+1))-(Y(i,j+1)-Y(i,j)));
    C(i-3,j-3)=1/hy*((Y(i,j+1)-Y(i,j))-(Y(i,j)-Y(i,j-1)));
    D(i-3,j-3)=1/hy*((Y(i,j)-Y(i,j-1))-(Y(i,j-1)-Y(i,j-2)));
    w(i-3,j-3)=1/(12*hy)*(-(Y(i,j-1)-Y(i,j-2))+7*(Y(i,j)-Y(i,j-1))+7*(Y(i,j+1)-Y(i,j))-(Y(i,j+2)-Y(i,j+1)))+weno(A(i-3,j-3),B(i-3,j-3),C(i-3,j-3),D(i-3,j-3));
end
end
end
end
