function [chi_x, chi_y] = cutoff_grad(epsilon,p,xs,ys)
x = p(:,1);  y = p(:,2);
r = sqrt((x-xs).^2 + (y-ys).^2);
a = epsilon;  b = 2*epsilon;
N = f(b-r); D = N+f(r-a);
dN = -f1(b-r); dD = dN+f1(r-a);
d = (dN.*D-N.*dD)./(D.*D);
drx = (x-xs)./r; dry = (y-ys)./r;
chi_x = d.*drx; chi_y = d.*dry;
chi_x(r<a) = 0; chi_y(r<a) = 0;


function f = f(x)
f = (x>0).*exp(-1./abs(x));

function f1 = f1(x)  % derivative of f(x)
f1 = (x>0).*f(x)./(x.*x);
