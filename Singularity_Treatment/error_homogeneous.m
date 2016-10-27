%% Compute errors of Ray-FEM solution in homogeneous medium
% 1. need suitable parameters: epsilon, omega, wpml
% 2. fix epsilon error \sim \omega^??

function error = error_homogeneous(omega,epsilon)

xs = 0;   ys = 0;                  % point source location
speed = @(x) ones(size(x,1),1);    % medium speed

fquadorder = 3;                    % numerical quadrature order

wavelength = 2*pi/omega; 
NPW = 6;
h = wavelength/NPW;
h = 1/round(1/h);
   
wpml = 0.08;%3*round(wavelength/h)*h;                 % width of PML
sigmaMax = 25/wpml;                % Maximun absorbtion

a = 1/2;
[node,elem] = squaremesh([-a,a,-a,a],h);
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
du(x>=a-wpml)=0; du(x<= -a+wpml) = 0;
du(y>=a-wpml)=0; du(y<= -a+wpml) = 0;
r = sqrt(x.*x+y.*y);
du(r<h/2) = 0;


%% Compute the errors of wavefield at y=0.4
px = 0.4;
px = round(px/h)*h;
indx = 1 + round((px+a)/h);
n = round(sqrt(N));
du = reshape(du,n,n);
uex = reshape(ub,n,n);

du = du(indx,:);
eu = uex(indx,:);

max_err = norm(du,inf);
rel_max_err = norm(du,inf)/norm(eu,inf);
l2_err = norm(du)*h;
rel_l2_err = norm(du)/norm(eu);
error = [max_err,rel_max_err,l2_err,rel_l2_err];





