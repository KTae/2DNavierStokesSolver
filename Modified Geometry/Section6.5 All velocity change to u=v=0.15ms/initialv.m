% build 2D initial v array
function v=initialv(x,y)
v=zeros(length(x),length(y));
for j=1:length(y)
    for i=1:length(x)
%         v(i,j)=cos(1.4*y(j))+cos(1.3*x(i));
        v(i,j)=0;
    end
end
end