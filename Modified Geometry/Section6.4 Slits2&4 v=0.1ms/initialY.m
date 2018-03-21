% build 2D initial Y array
function Y=initialY(x,y)
Y=zeros(length(x),length(y));
for j=1:length(y)
    for i=1:length(x)
%         Y(i,j)=sin(0.7*x(i))+cos(1.2*y(j));
              Y(i,j)=0;
    end
end
end