% project/correct velocities using Lagrange multiplier
function [unew,vnew]=projection(u,v,phi,dt,dx,dy,yu,xv,m,n)
unew=zeros(m+1,n+2);
vnew=zeros(m+2,n+1);
for j=2:n+1
    for i=2:m
unew(i,j)=u(i,j)-dt/dx*(phi(i+1,j)-phi(i,j));
    end
end
for j=2:n
    for i=2:m+1
   vnew(i,j)=v(i,j)-dt/dy*(phi(i,j+1)-phi(i,j));
    end
end
% apply boundary
unew=ubc(unew,yu);
vnew=vbc(vnew,xv);
end