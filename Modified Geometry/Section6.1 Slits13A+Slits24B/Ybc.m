% apply bc to Y input 106x26 output 106x26 
function Y=Ybc(Y,m,n,x,y)
% outlet bc
Y(m+4,:)=Y(m+3,:);
Y(m+5,:)=Y(m+2,:);
Y(m+6,:)=Y(m+1,:);
% inlet bc
for i=1:length(y)
if y(i)<0.5
    Y(3,i)=2-Y(4,i);
    Y(2,i)=2-Y(5,i);
    Y(1,i)=2-Y(6,i);
else 
    Y(3,i)=-Y(4,i);
    Y(2,i)=-Y(5,i);
    Y(1,i)=-Y(6,i);
end
end
% wall bc
Y(:,3)=Y(:,4);
Y(:,2)=Y(:,5);
Y(:,1)=Y(:,6);
Y(:,n+4)=Y(:,n+3);
Y(:,n+5)=Y(:,n+2);
Y(:,n+6)=Y(:,n+1);
% slits
for i=1:length(x)
    if x(i)<=1.125&&x(i)>=0.875
        Y(i,3)=2-Y(i,4);
        Y(i,2)=2-Y(i,5);
        Y(i,1)=2-Y(i,6);
    else if x(i)<=3.125&&x(i)>=2.875
         Y(i,3)=2-Y(i,4);
         Y(i,2)=2-Y(i,5);
         Y(i,1)=2-Y(i,6);
        else if x(i)>=1.875&&x(i)<=2.125
                Y(i,n+6)=-Y(i,n+1);
                Y(i,n+5)=-Y(i,n+2);
                Y(i,n+4)=-Y(i,n+3);
            else if x(i)>=3.875&&x(i)<=4.125         
                Y(i,n+6)=-Y(i,n+1);
                Y(i,n+5)=-Y(i,n+2);
                Y(i,n+4)=-Y(i,n+3);
                end
            end
        end
    end
end 
end