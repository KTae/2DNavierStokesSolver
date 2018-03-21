%calculate the rhs for possion equation f(m+2,n+2)
function [f]=possionrhs(u,v,dx,dy,dt,m,n)
f=zeros(m+2,n+2);
for j=2:n+1
    for i=2:m+1
   f(i,j)=1/dt*((u(i,j)-u(i-1,j))/dx+(v(i,j)-v(i,j-1))/dy);
    end
end
end