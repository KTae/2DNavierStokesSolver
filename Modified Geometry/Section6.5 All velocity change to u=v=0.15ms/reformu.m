% move u to the cell centers
function unew=reformu(u,m,n)
unew=zeros(m,n);
for j=1:n
    for i=1:m
        unew(i,j)=0.5*(u(i,j+1)+u(i+1,j+1));
    end
end
end