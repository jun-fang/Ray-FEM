%% Compute errors of Ray-FEM solution in homogeneous medium
% 1. need suitable parameters: epsilon, omega, wpml
% 2. fix epsilon error \sim \omega^??



% need to figure out the reason of convergence rate
% 1. check the source point
% 2. epsilon
% 3. error region
% 4. wpml
function error = error_homogeneous(omega,epsilon)

xs = 0.1;   ys = 0.1;                  % point source location
% xs = 0;   ys = 0;
speed = @(x) ones(size(x,1),1);    % medium speed

fquadorder = 3;                    % numerical quadrature order

wavelength = 2*pi/omega; 
NPW = 6;
h = wavelength/NPW;
h = 1/round(1/h);
   
wpml = 0.1;%3*round(wavelength/h)*h;                 % width of PML
sigmaMax = 25/wpml;                % Maximun absorbtion

a = 0.5;
% [node,elem] = squaremesh([-a,a,-a,a],h);
[node,elem] = squaremesh([-.4,.6,-.4,.6],h);
Nray = 1; N = size(node,1); Ndof = N;

%% Exact ray information
xx = node(:,1)-xs;  yy = node(:,2)-ys;
rr = sqrt(xx.^2 + yy.^2);
ray = atan2(yy,xx);
ray = exp(1i*ray).*(rr>10*eps);

tic;
%% Assembling
[A] = assemble_Helmholtz_matrix_with_ray_1(node,elem,omega,wpml,sigmaMax,speed,ray,fquadorder);
b = assemble_RHS_with_ray_sing(node,elem,epsilon,omega,speed,xs,ys,ray,fquadorder);


%% Boundaries
[~,~,isBdNode] = findboundary(elem);
rep_isBdNode = repmat(isBdNode,1,Nray);
isBdNode = rep_isBdNode(:);
freeNode = find(~isBdNode);


%% Solve the linear system Au = b
v = zeros(Ndof,1);
v(freeNode) = A(freeNode,freeNode)\b(freeNode);


%% Compute solution values at grid nodes
grad = ray(:);
grad = [real(grad),imag(grad)];
repnode = repmat(node,Nray,1);
temp = grad(:,1).*repnode(:,1) + grad(:,2).*repnode(:,2);

k = omega./speed(node);           % wavenumber
kk = repmat(k,1,Nray);
u = v.*exp(1i*kk(:).*temp);
u = reshape(u,N,Nray);
u = sum(u,2);


xx = node(:,1)-xs;  yy = node(:,2)-ys;
rr = sqrt(xx.^2 + yy.^2);
ub = 1i/4*besselh(0,1,omega*rr);
aa = epsilon; bb = 2*epsilon;
x_eps = cutoff(aa,bb,node,xs,ys);
v = (1-x_eps).*ub;
v(rr<aa)=0;

du = u-v;
x = node(:,1); y = node(:,2);
du(x>=max(x)-wpml)=0; du(x<= min(x)+wpml) = 0;
du(y>=max(y)-wpml)=0; du(y<= min(y)+wpml) = 0;

max_err = norm(du,inf);
rel_max_err = norm(du,inf)/norm(v,inf);
l2_err = norm(du)*h;
rel_l2_err = norm(du)/norm(v);

rhs = sing_rhs(epsilon,omega,node,xs,ys);
rhsL2 = norm(rhs)*h;
error = [max_err,rel_max_err,l2_err,rel_l2_err,rhsL2];

if(0)


%% Construct the solution
xx = node(:,1)-xs;  yy = node(:,2)-ys;
rr = sqrt(xx.^2 + yy.^2);
ub = 1i/4*besselh(0,1,omega*rr);
aa = epsilon; bb = 2*epsilon;
x_eps = cutoff(aa,bb,node,xs,ys);
v = ub.*x_eps;
us = u;
uu = us + v;
toc;


du = ub -uu;
x = node(:,1); y = node(:,2);
xx = x -xs; yy = y - ys;
du(x>=max(x)-wpml)=0; du(x<= min(x)+wpml) = 0;
du(y>=max(y)-wpml)=0; du(y<= min(y)+wpml) = 0;
r = sqrt(xx.*xx+yy.*yy);
du(r<h/2) = 0;


%% Compute the errors of wavefield at y=0.4
% px = 0.4;
% px = round(px/h)*h;
% indx = 1 + round((px+a)/h);
% n = round(sqrt(N));
% du = reshape(du,n,n);
% uex = reshape(ub,n,n);
% 
% du = du(indx,:);
% eu = uex(indx,:);


%% Compute the errors in a box [0.2,0.4]^2
du = ub -uu;
x = node(:,1); y = node(:,2);
du = du((x>=0.2-eps)&(x<=0.4+eps)&(y>=0.2-eps)&(y<=0.4+eps));
eu = ub((x>=0.2-eps)&(x<=0.4+eps)&(y>=0.2-eps)&(y<=0.4+eps));


%% Compute the errors in annulus epsilon < r < 2*epsilon
% du = du(r>=aa & r<=bb);
% eu = ub(r>=aa & r<=bb);


%% Compute the errors outside one wavelength to the source point
% idx_pml = (abs(x)<=max(x)-wpml)&(abs(y)<=max(y)-wpml);
% du = du(idx_pml & (r>wavelength));
% eu = ub(idx_pml & (r>wavelength));

max_err = norm(du,inf);
rel_max_err = norm(du,inf)/norm(eu,inf);
l2_err = norm(du)*h;
rel_l2_err = norm(du)/norm(eu);


rhs = sing_rhs(epsilon,omega,node,xs,ys);
rhsL2 = norm(rhs)*h;
error = [max_err,rel_max_err,l2_err,rel_l2_err,rhsL2];

%figure()

end

