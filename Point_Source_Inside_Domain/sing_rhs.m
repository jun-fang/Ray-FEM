function rhs = sing_rhs(epsilon,omega,p)
r = sqrt(p(:,1).^2 + p(:,2).^2);
a = 1/epsilon;  b = 2/epsilon;
ub = -1i/4*besselh(0,1,omega*r);
ub_g1 = 1i/4*omega*besselh(1,1,omega*r)./r.*p(:,1);
ub_g2 = 1i/4*omega*besselh(1,1,omega*r)./r.*p(:,2);
% ub_grad = [ub_g1, ub_g2];

cf_grad = cutoff_gradient(epsilon,p);
cf_lap = cutoff_laplacian(epsilon,p);

rhs = -2*(ub_g1.*cf_grad(:,1) + ub_g2.*cf_grad(:,2)) - ub.*cf_lap;
rhs(r<=a) = 0;
rhs(r>=b) = 0;