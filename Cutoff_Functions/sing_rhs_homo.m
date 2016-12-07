%% Singularity treatment right hand side function 
% homogeneous medium in the support of cut-off function

function rhs = sing_rhs_homo(epsilon,omega,p,xs,ys)

a = epsilon;  b = 2*epsilon;
x = (p(:,1)-xs);  y = (p(:,2)-ys);
r = sqrt(x.^2 + y.^2);

ub = 1i/4*besselh(0,1,omega*r);                  % Babich expression
ub_g1 = -1i/4*omega*besselh(1,1,omega*r)./r.*x;  % partial derivative wrt x
ub_g2 = -1i/4*omega*besselh(1,1,omega*r)./r.*y;  % partial derivative wrt y

cf_grad = cutoff_gradient(a,b,p,xs,ys);
cf_lap = cutoff_laplacian(a,b,p,xs,ys);

rhs = 2*(ub_g1.*cf_grad(:,1) + ub_g2.*cf_grad(:,2)) + ub.*cf_lap;
rhs(r<a) = 0;   rhs(r>b) = 0;