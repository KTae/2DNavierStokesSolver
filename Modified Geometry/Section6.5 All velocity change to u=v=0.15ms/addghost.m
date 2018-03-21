% this function is to add ghost cell values to T
function [T]=addghost(T)
T(1,:)=T(2,:);
T(end,:)=T(end-1,:);
T(:,1)=T(:,2);
T(:,end)=T(:,end-1);
end