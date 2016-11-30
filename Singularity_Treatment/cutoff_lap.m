function cfl = cutoff_lap(epsilon,p,xs,ys)
x = p(:,1);  y = p(:,2);
r = sqrt((x-xs).^2 + (y-ys).^2);
a = epsilon;  b = 2*epsilon;
N = f(b-r); D = N+f(r-a);
gr2 = ones(size(r)); lr = 1./r;
gDgN = f1(b-r).*( f1(b-r)-f1(r-a) ).*gr2;
gDgD = ( f1(b-r)-f1(r-a) ).^2.*gr2;
lD = ( f2(b-r)+f2(r-a) ).*gr2 + ( f1(r-a)-f1(b-r) ).*lr;
lN = f2(b-r).*gr2 - f1(b-r).*lr;

numerator = (D.*lN-N.*lD).*D - 2*D.*gDgN + 2*N.*gDgD;
denominator = D.^3;
cfl = numerator./denominator;
cfl(r<a+10*eps)=0;

function f = f(x)
f = (x>0).*exp(-1./abs(x));

function f1 = f1(x)  % derivative of f(x)
f1 = (x>0).*f(x)./(x.*x);

function f2 = f2(x)
f2 = (x>0).*f(x).*(x.^(-4) - 2*x.^(-3));
