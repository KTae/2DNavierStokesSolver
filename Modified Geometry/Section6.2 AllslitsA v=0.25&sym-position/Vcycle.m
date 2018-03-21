function [phi,err]=Vcycle(phi,f,m,h,options,nmax,alpha)
% set initial parameters
p=log2(m/25)+1;
% judge m
if rem(p,1)~=0
    err=2;
    warning('the number of elements is not a correct value')
    return
end
% err is 1 or 0 need to determine later
if rem(p,1)==0
% err=1;
% err=0;
    
M=zeros(1,p); 
M(1)=m;
% build Mesh space array M
for i=2:p
M(i)=M(i-1)/2;
end
% build N array
N=M/5;
% build h array
H=zeros(1,p);
H(1)=h;
for i=2:p
H(i)=2*H(i-1);
end


% option 1
if options==1
% initial processing
r=zeros(M(1)+2,N(1)+2,length(M));
phi(:,:,1)=phi;
% loop
for j=1:nmax
eps=zeros(M(1)+2,N(1)+2,length(M));
phi(:,:,1)=GaussSeidel(phi(:,:,1),f,H(1),M(1),N(1));
r(:,:,1)=calcResidual(phi(:,:,1),f,H(1),M(1),N(1));
for q=2:p
    rhs(1:M(q)+2,1:N(q)+2,q)=restrict(r(:,:,q-1),M(q),N(q));
    eps(:,:,q)=GaussSeidel(eps(:,:,q),rhs(:,:,q),H(q),M(q),N(q));
    r(1:M(q)+2,1:N(q)+2,q)=calcResidual(eps(:,:,q),rhs(:,:,q),H(q),M(q),N(q));
end
epsc=zeros(M(1)+2,N(1)+2,p);
for q=p-1:-1:2
    epsc(1:M(q)+2,1:N(q)+2,q)=prolong(eps(1:M(q+1)+2,1:N(q+1)+2,q+1),M(q+1),N(q+1));
    eps(:,:,q)=correct(eps(:,:,q),epsc(:,:,q));
    eps(:,:,q)=GaussSeidel(eps(:,:,q),rhs(:,:,q),H(q),M(q),N(q));
end
epsc(:,:,1)=prolong(eps(:,:,2),M(2),N(2));
phi(:,:,1)=correct(phi(:,:,1),epsc(:,:,1));
end
err=0;
end

% option 2
if options==2
% initial processing
r=zeros(M(1)+2,N(1)+2,length(M));
error=zeros(M(1)+2,N(1)+2);
Linf=zeros(nmax,1);
phi(:,:,1)=phi;
% loop
for j=1:nmax
eps=zeros(M(1)+2,N(1)+2,length(M));
phi(:,:,1)=GaussSeidel(phi(:,:,1),f,H(1),M(1),N(1));
r(:,:,1)=calcResidual(phi(:,:,1),f,H(1),M(1),N(1));
for q=2:p
    rhs(1:M(q)+2,1:N(q)+2,q)=restrict(r(:,:,q-1),M(q),N(q));
    eps(:,:,q)=GaussSeidel(eps(:,:,q),rhs(:,:,q),H(q),M(q),N(q));
    r(1:M(q)+2,1:N(q)+2,q)=calcResidual(eps(:,:,q),rhs(:,:,q),H(q),M(q),N(q));
end
epsc=zeros(M(1)+2,N(1)+2,p);
for q=p-1:-1:2
    epsc(1:M(q)+2,1:N(q)+2,q)=prolong(eps(1:M(q+1)+2,1:N(q+1)+2,q+1),M(q+1),N(q+1));
    eps(:,:,q)=correct(eps(:,:,q),epsc(:,:,q));
    eps(:,:,q)=GaussSeidel(eps(:,:,q),rhs(:,:,q),H(q),M(q),N(q));
end
epsc(:,:,1)=prolong(eps(:,:,2),M(2),N(2));
phi(:,:,1)=correct(phi(:,:,1),epsc(:,:,1));
% calculate error
for jj=2:length(phi(1,:))-1
for i=2:length(phi(:,1))-1
error(i,jj)=f(i,jj)-1/H(1)^2*(phi(i-1,jj)+phi(i+1,jj)+phi(i,jj+1)+phi(i,jj-1)-4*phi(i,jj));
end
end
% calculate infinity norm of error
Linf(j)=max(max(abs(error)));
% absolute convergence determine below max n
if Linf(j)<alpha
err=0;
return
end
end
% absolute convergence determine after max n
if Linf(end)>alpha
    err=1;
    warning ('set iteration times did not satisfying criterium 2');
end
if Linf(end)==alpha
    err=1;
    warning ('set iteration times did not satisfying criterium 2');
end
end


% option 3
if options==3
% initial processing
r=zeros(M(1)+2,N(1)+2,length(M));
error=zeros(M(1)+2,N(1)+2);
Linf=zeros(nmax,1);

% calculate infinity norm of initial guess Linf0
residual0=zeros(M(1)+2,N(1)+2);
for i=2:length(phi(:,1))-1
for j=2:length(phi(1,:))-1
residual0(i,j)=f(i,j)-1/H(1)^2*(phi(i-1,j)+phi(i,j+1)+phi(i+1,j)+phi(i,j-1)-4*phi(i,j));
end
end
Linf0=max(max(abs(residual0)));

phi(:,:,1)=phi;
% loop
for j=1:nmax
eps=zeros(M(1)+2,N(1)+2,length(M));
phi(:,:,1)=GaussSeidel(phi(:,:,1),f,H(1),M(1),N(1));
r(:,:,1)=calcResidual(phi(:,:,1),f,H(1),M(1),N(1));
for q=2:p
    rhs(1:M(q)+2,1:N(q)+2,q)=restrict(r(:,:,q-1),M(q),N(q));
    eps(:,:,q)=GaussSeidel(eps(:,:,q),rhs(:,:,q),H(q),M(q),N(q));
    r(1:M(q)+2,1:N(q)+2,q)=calcResidual(eps(:,:,q),rhs(:,:,q),H(q),M(q),N(q));
end
epsc=zeros(M(1)+2,N(1)+2,p);
for q=p-1:-1:2
    epsc(1:M(q)+2,1:N(q)+2,q)=prolong(eps(1:M(q+1)+2,1:N(q+1)+2,q+1),M(q+1),N(q+1));
    eps(:,:,q)=correct(eps(:,:,q),epsc(:,:,q));
    eps(:,:,q)=GaussSeidel(eps(:,:,q),rhs(:,:,q),H(q),M(q),N(q));
end
epsc(:,:,1)=prolong(eps(:,:,2),M(2),N(2));
phi(:,:,1)=correct(phi(:,:,1),epsc(:,:,1));
% calculate error
for jj=2:length(phi(1,:))-1
for i=2:length(phi(:,1))-1
error(i,jj)=f(i,jj)-1/H(1)^2*(phi(i-1,jj)+phi(i+1,jj)+phi(i,jj+1)+phi(i,jj-1)-4*phi(i,jj));
end
end
% calculate infinity norm of error
Linf(j)=max(max(abs(error)));
% calculate the ratio
ratio(j)=Linf(j)/Linf0;
% relative convergence in the middle iterations
if ratio(j)<alpha
    err=0;
    return
end
end
% judge the last time iteration
if ratio(end)>alpha
    err=1;
   warning ('set iteration times did not satisfying criterium 3');
end
if ratio(end)==alpha
    err=1;
   warning ('set iteration times did not satisfying criterium 3');
end
end

end
end







