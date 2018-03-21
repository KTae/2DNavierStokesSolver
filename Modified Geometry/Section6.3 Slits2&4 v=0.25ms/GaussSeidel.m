function y=GaussSeidel(phi,rhs,h,m,n)
% GS method 
for j=2:m+1
    for i=2:n+1
     phi(j,i)=0.25*(phi(j-1,i)+phi(j+1,i)+phi(j,i+1)+phi(j,i-1))-0.25*h^2*rhs(j,i);
    end  
end
% apply ghost cell values
phi(1,:)=phi(2,:);
phi(m+2,:)=phi(m+1,:);
phi(:,1)=phi(:,2);
phi(:,n+2)=phi(:,n+1);
y=phi;


