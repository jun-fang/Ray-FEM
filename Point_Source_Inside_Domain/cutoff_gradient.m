function grad = cutoff_gradient(epsilon,p)
x = p(:,1);  y = p(:,2);
r = sqrt(x.^2 + y.^2);
a = 1/epsilon;  b = 2/epsilon;
temp = f(b-r) + f(r-a);
g1 = f1(b-r).*(-x./r).*temp - f(b-r).*( f1(b-r).*(-x./r) + f1(r-a).*x./r );
g2 = f1(b-r).*(-y./r).*temp - f(b-r).*( f1(b-r).*(-y./r) + f1(r-a).*y./r );
g1 = g1./(temp.*temp);
g2 = g2./(temp.*temp);
grad = [g1, g2];

function f = f(x)
f = (x>10*eps).*exp(-1./x);

function f1 = f1(x)
f1 = f(x)./(x.*x);
