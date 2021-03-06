% GCI analysis
R=[0.0516978946635019;0.0530658609520460;0.0571541265276045];
% GCI analysis
% initial processing
Fsec=1.25;
r=2;
% calculate order of convergence p
p=log(((R(3)-R(2))/(R(2)-R(1))))/log(r);
% determine estimate for exact solution using Richardson extrapolation
f0=R(1)+(R(1)-R(2))/(r^p-1);
% Grid Convergence Index
e12=abs((R(1)-R(2))/R(1));
GCI21=Fsec*e12/(r^p-1);
e23=abs((R(2)-R(3))/R(2));
GCI23=Fsec*e23/(r^p-1);
% Asymptotic range of convergence
k=GCI21/GCI23*r^p;
% final answer output in command window
per=abs(GCI21*100);
formatSpec = 'The value of R at t=0.03s is %f0 +/- %G%%';
GCI_analysis_result = sprintf(formatSpec,f0,per)