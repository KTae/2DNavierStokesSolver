function y=calcResidual(phi,rhs,h,m,n)
% calculate residual
y=zeros(m+2,n+2);
for j=2:n+1
for i=2:m+1
    y(i,j)=rhs(i,j)-1/h^2*(phi(i+1,j)+phi(i-1,j)+phi(i,j-1)+phi(i,j+1)-4*phi(i,j));
end
end
% set ghost cell values
y(1,:)=y(2,:);
y(end,:)=y(end-1,:);
y(:,1)=y(:,2);
y(:,end)=y(:,end-1);
end