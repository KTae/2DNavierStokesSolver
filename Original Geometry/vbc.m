% add bc to v input(m+2,n+1) output(m+2,n+1)
function m=vbc(m,x)
% upper and lower wall bc
m(:,1)=vslitlow(x);
m(:,end)=vslitup(x);
% inlet bc (Drichilet bc)
m(1,:)=-m(2,:);
% outlet bc (Drichilet bc)
m(end,:)=-m(end-1,:);
end

%  add v velocity to lower slits
function m=vslitlow(x)
m=zeros(1,length(x));
for i=1:length(x)
    if x(i)>=0.875&&x(i)<=1.125
        m(i)=0.5;
    else if x(i)>=2.875&&x(i)<=3.125
        m(i)=0.5;
        end
    end
end
end

%  add v velocity to upper slits
function m=vslitup(x)
m=zeros(1,length(x));
for i=1:length(x)
    if x(i)>=1.875&&x(i)<=2.125
        m(i)=-0.5;
    else if x(i)>=3.875&&x(i)<=4.125
        m(i)=-0.5;
        end
    end
end
end