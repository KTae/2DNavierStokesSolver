% add bc to u  input(m+1,n+2) output(m+1,n+2)
function m=ubc(m,y)
% left bc (inlet source A,B) 
m(1,:)=uinlet(y);
% right bc(outlet source 0st order in dubug)
m(end,:)=m(end-1,:);
% upper and lower wall bc
m(:,1)=-m(:,2);
m(:,end)=-m(:,end-1);
end

% inlet of u for initial condition
function uy=uinlet(y)
uy=zeros(1,length(y));
for i=1:length(y)
if y(i)<=0.5
    uy(i)=-16*(y(i)-0.25)^2+1;
else
    uy(i)=-16*(y(i)-0.75)^2+1;
end
end
end