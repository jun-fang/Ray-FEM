%% gradient of the cut-off function
% output is 0 if r>=b or r<=a, smooth in the intermediate annulus

function grad = cutoff_gradient(a,b,p,xs,ys)
x = p(:,1)-xs;  y = p(:,2)-ys;
r = sqrt(x.^2 + y.^2);
temp = f(b-r) + f(r-a);
g1 = f1(b-r).*(-x./r).*temp - f(b-r).*( f1(b-r).*(-x./r) + f1(r-a).*x./r );
g2 = f1(b-r).*(-y./r).*temp - f(b-r).*( f1(b-r).*(-y./r) + f1(r-a).*y./r );
g1 = g1./(temp.*temp);  g2 = g2./(temp.*temp);
g1(r<a+eps) = 0;  g1(r>b-eps) = 0;
g2(r<a+eps) = 0;  g2(r>b-eps) = 0;
grad = [g1, g2];


function f = f(x)
f = (x>0).*exp(-1./abs(x));

function f1 = f1(x)  % derivative of f(x)
f1 = f(x)./(x.*x);
