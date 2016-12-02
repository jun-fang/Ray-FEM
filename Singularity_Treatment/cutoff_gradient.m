%% gradient of the cut-off function
% output is 0 if r>=b or r<=a, smooth in the annulus {r| a<r<b}
% optional choice for a,b: a = epsilon;  b = 2*epsilon;

function cfg = cutoff_gradient(a,b,p,xs,ys)
x = p(:,1);  y = p(:,2);
r = sqrt((x-xs).^2 + (y-ys).^2);    % distance to (xs,ys)
N = f(b-r); D = N+f(r-a);           % numerator and denominator of cut-off function
dN = -f1(b-r); dD = dN+f1(r-a);     % differentiate numerator and denominator
d = (dN.*D-N.*dD)./(D.*D);          % quotient rule
drx = (x-xs)./r; dry = (y-ys)./r;   % partial derivative wrt x and y
chi_x = d.*drx; chi_y = d.*dry;     % gradient
chi_x(r<a) = 0; chi_y(r<a) = 0;     % eliminate Nan at (xs,ys)
cfg = [chi_x, chi_y];

function f = f(x)
f = (x>0).*exp(-1./abs(x));

function f1 = f1(x)  % derivative of f(x)
f1 = (x>0).*f(x)./(x.*x);
