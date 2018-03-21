function rcoarse=restrict(rfine,mcoarse,ncoarse)
% restrict the residual to coarse mesh 
rcoarse=zeros(mcoarse+2,ncoarse+2);
for j=2:ncoarse+1
for i=2:mcoarse+1
    rcoarse(i,j)=0.25*(rfine(2*i-2,2*j-2)+rfine(2*i-2,2*j-1)+rfine(2*i-1,2*j-2)+rfine(2*i-1,2*j-1));
end
end
% add ghost cell values
rcoarse(1,:)=rcoarse(2,:);
rcoarse(end,:)=rcoarse(end-1,:);
rcoarse(:,1)=rcoarse(:,2);
rcoarse(:,end)=rcoarse(:,end-1);
end