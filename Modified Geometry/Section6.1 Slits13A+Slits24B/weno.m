% calculate the coefficient in WENO-5
function phi=weno(a,b,c,d)
e=1e-6;
IS0=13*(a-b)^2+3*(a-3*b)^2;
IS1=13*(b-c)^2+3*(b+c)^2;
IS2=13*(c-d)^2+3*(3*c-d)^2;
a0=1/(e+IS0)^2;
a1=6/(e+IS1)^2;
a2=3/(e+IS2)^2;
w0=a0/(a0+a1+a2);
w2=a2/(a0+a1+a2);
phi=1/3*w0*(a-2*b+c)+1/6*(w2-0.5)*(b-2*c+d);
end