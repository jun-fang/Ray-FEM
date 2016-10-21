function x_eps = cutoff(epsilon,p)
r = sqrt(p(:,1).^2 + p(:,2).^2);
a = 1/epsilon;
b = 2/epsilon;
x_eps = g(r,a,b);

% function f = f(p)
% r = p(:,1).^2 + p(:,2).^2;
% f = (r<1).*exp(-1./(1-r));
% 
% function phi = phi(p)
% phi = (f(p)>0).*f(p)./(f(p) + f(2-p));

function f = f(x)
f = (x>10*eps).*exp(-1./x);

function g = g(x,a,b)
g = f(b-x)./( f(b-x)+f(x-a) );

