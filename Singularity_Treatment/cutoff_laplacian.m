%% Laplacian of the cut-off function
% output is 0 if r>=b or r<=a, smooth in the intermediate annulus

function cfl = cutoff_laplacian(a,b,p,xs,ys)
x = p(:,1)-xs;  y = p(:,2)-ys;
r = sqrt(x.^2 + y.^2);
cfl = g2(r,a,b) + g1(r,a,b)./r;
cfl(r<a+eps) = 0;   cfl(r>b-eps) = 0;

function f = f(x)
f = (x>0).*exp(-1./abs(x));

function f1 = f1(x)
f1 = f(x)./(x.*x);

function f2 = f2(x)
f2 = f(x).*(x.^(-4) - 2*x.^(-3));

function g = g(x,a,b)
g = f(b-x)./( f(b-x)+f(x-a) );

function g1 = g1(x,a,b)
g1 = -( f1(b-x).*(f(b-x) + f(x-a)) + f(b-x).*( f1(x-a) - f1(b-x) ) )...
    ./( f(b-x) + f(x-a) ).^2;

function g2 = g2(x,a,b)
temp = f(b-x)+f(x-a);
temp0 = f1(b-x).*(f(b-x)+f(x-a)) + f(b-x).*(f1(x-a)-f1(b-x));
temp1 = f2(b-x).*f(x-a) - f(b-x).*f2(x-a);
g2 = temp1./(temp.^2) + 2*temp0./(temp.^3).*(f1(x-a) - f1(b-x));
