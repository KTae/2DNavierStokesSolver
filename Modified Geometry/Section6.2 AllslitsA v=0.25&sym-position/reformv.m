% move v to the cell centers
function vnew=reformv(v,m,n)
vnew=zeros(m,n);
for j=1:n
    for i=1:m
        vnew(i,j)=0.5*(v(i+1,j)+v(i+1,j+1));
    end
end
end