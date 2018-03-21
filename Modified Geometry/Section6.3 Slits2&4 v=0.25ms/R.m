% calculate R(t)
function M=R(Y,m,n,dx,dy)
r=0;
for  j=4:n-3
    for i=4:m-3
        r=r+Y(i,j)*(1-Y(i,j))*dx*dy;
    end
end
M=r/5;