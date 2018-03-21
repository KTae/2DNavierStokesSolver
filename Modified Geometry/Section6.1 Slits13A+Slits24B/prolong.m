function rfine=prolong(rcoarse,mcoarse,ncoarse)
% prolong the residual on coarse mesh
rfine=zeros(2*mcoarse+2,2*ncoarse+2);
for j=2:ncoarse+1
for i=2:mcoarse+1
    rfine(2*i-1,2*j-1)=rcoarse(i,j);
    rfine(2*i-1,2*j-2)=rcoarse(i,j);
    rfine(2*i-2,2*j-1)=rcoarse(i,j);
    rfine(2*i-2,2*j-2)=rcoarse(i,j);
end
end
% add ghost cell values
rfine(1,:)=rfine(2,:);
rfine(end,:)=rfine(end-1,:);
rfine(:,1)=rfine(:,2);
rfine(:,end)=rfine(:,end-1);
end