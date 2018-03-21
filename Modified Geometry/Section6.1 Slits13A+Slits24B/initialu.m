% build 2D initial u array
function u=initialu(x,y)
u=zeros(length(x),length(y));
for j=1:length(y)
    for i=1:length(x)
%         u(i,j)=sin(1.2*x(i)+0.1)+sin(1.4*y(j));
        u(i,j)=0;
    end
end
end