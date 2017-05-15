function rhs = singularity_RHS(epsilon,xs,ys,p,ub,ub_g1,ub_g2)
%% Singularity treatment right hand side function 
%  
%         f = 2 \nabla u_b \dot \nabla \chi_{\epsilon} + u_b \Delta \chi_{\epsilon}
%                
%  INPUT: 
%      epsilon: cut-off supoort size 
%      (xs, ys): source location
%      p: nx2 arrray, x, y coordinates of n nodes
%      ub: nx1 vector, Babich expression
%      ub_g1: nx1 vector, partial derivative \patial ub / \partial x
%      ub_g2: nx1 vector, partial derivative \patial ub / \partial y
%
%   OUTPUT:
%      rhs = f: nx1 vector
% 

a = epsilon;  b = 2*epsilon;
cf_grad = cutoff_gradient(a,b,p,xs,ys);
cf_lap = cutoff_laplacian(a,b,p,xs,ys);

rhs = 2*(ub_g1(:).*cf_grad(:,1) + ub_g2(:).*cf_grad(:,2)) + ub(:).*cf_lap;
x = (p(:,1)-xs);  y = (p(:,2)-ys);
r = sqrt(x.^2 + y.^2);
rhs(r<a) = 0;   rhs(r>b) = 0;


%% homogeneous medium case
% x = (p(:,1)-xs);  y = (p(:,2)-ys);
% r = sqrt(x.^2 + y.^2);
% ub = 1i/4*besselh(0,1,omega*r);                  % Babich expression
% ub_g1 = -1i/4*omega*besselh(1,1,omega*r)./r.*x;  % partial derivative wrt x
% ub_g2 = -1i/4*omega*besselh(1,1,omega*r)./r.*y;  % partial derivative wrt y

