%% Cut-off function
% output is 0 if r>=b and 1 if r<=a, smooth in the intermediate annulus

function cf = cutoff(a,b,p,xs,ys)
r = sqrt((p(:,1)-xs).^2 + (p(:,2)-ys).^2);
cf = g(r,a,b);
cf(r<a+eps) = 1;   cf(r>b-eps) = 0;


function f = f(x)
f = (x>0).*exp(-1./abs(x));   
% It is important to use absolute value, otherwise it may produce 'Nan' 
% for negative numbers close to 0, like -0.001

function g = g(x,a,b)
g = f(b-x)./( f(b-x)+f(x-a) );

