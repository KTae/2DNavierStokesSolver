function d = GE(a,b,c,d)
% This function is to solve the tri-diagonal matrix using Gauss
% eliminatioon method. 
% The inputs are 4 arrays, a,b,c,d which can stand for the whole
% tri-diagonal matrix.
% b is the vector of diagonal entries, a is the array of numbers below
% diagonal entries and c is the numbers above the diagonal entries.
% The output is the vector of the correct solution of this given equation 
% restored in array d.
% However,this method will destory the original data vector b and d.
format long;
% In order to have more precision.
n=length(b);
% This is just the exactly derivation of Gauss elimination.
% Step 1: elimination
for i=2:n
    b(i)=b(i)-c(i-1)*a(i)/b(i-1);
    d(i)=d(i)-d(i-1)*a(i)/b(i-1);
end
% Step 2: back substitution
d(n)=d(n)/b(n);
for j=n-1:-1:1
    d(j)=(d(j)-c(j)*d(j+1))/b(j);
end
end

