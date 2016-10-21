% function x_eps = cutoff(epsilon,p)
% r = sqrt(p(:,1).^2 + p(:,2).^2);
% a = 1/epsilon;
% b = 2/epsilon;
% x_eps = g(r,a,b);

function cf = cutoff(a,b,p)
r = sqrt(p(:,1).^2 + p(:,2).^2);
cf = g(r,a,b);
cf(r<a+eps) = 1;   cf(r>b-eps) = 0;


function f = f(x)
f = (x>0).*exp(-1./x);

function g = g(x,a,b)
g = f(b-x)./( f(b-x)+f(x-a) );

