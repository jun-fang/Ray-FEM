%% Laplacian of the cut-off function
% output is 0 if r>=b or r<=a, smooth in the annulus {r| a<r<b}
% optional choice for a,b: a = epsilon;  b = 2*epsilon;

function cfl = cutoff_laplacian(a,b,p,xs,ys)

x = p(:,1);  y = p(:,2);
r = sqrt((x-xs).^2 + (y-ys).^2);     % distance to (xs,ys)
N = f(b-r); D = N+f(r-a);            % numerator and denominator of cut-off function
gr2 = ones(size(r)); lr = 1./r;      % |\nabla r| and \Delta r

gDgN = f1(b-r).*( f1(b-r)-f1(r-a) ).*gr2;     % \nabla D \dot \nabla N
gDgD = ( f1(b-r)-f1(r-a) ).^2.*gr2;           % \nabla D \dot \nabla D
lD = ( f2(b-r)+f2(r-a) ).*gr2 + ( f1(r-a)-f1(b-r) ).*lr;    % \Delta D
lN = f2(b-r).*gr2 - f1(b-r).*lr;                            % \Delta N

numerator = (D.*lN-N.*lD).*D - 2*D.*gDgN + 2*N.*gDgD;
denominator = D.^3;
cfl = numerator./denominator;
cfl(r<a+10*eps)=0;

function f = f(x)
f = (x>0).*exp(-1./abs(x));

function f1 = f1(x)  % first derivative of f(x)
f1 = (x>0).*f(x)./(x.*x);

function f2 = f2(x)  % second derivative of f(x)
f2 = (x>0).*f(x).*(x.^(-4) - 2*x.^(-3));
