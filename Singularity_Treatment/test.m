% xs = 0.1;   ys = 0.1;                  % point source location
xs = 0;   ys = 0;
speed = @(x) ones(size(x,1),1);    % medium speed
source = @(x) 1i*10*pi*exp(1i*10*pi*sqrt(x(:,1).^2+x(:,2).^2))./sqrt(x(:,1).^2+x(:,2).^2);
exu = @(x) exp(1i*10*pi*sqrt(x(:,1).^2+x(:,2).^2));

fquadorder = 3;                    % numerical quadrature order


omega = 10*pi;
epsilon = 0.131;
wavelength = 2*pi/omega; 

NPW = 12;
h = wavelength/NPW;
h = 1/round(1/h);
   
wpml = 0.1;%3*round(wavelength/h)*h;                 % width of PML
sigmaMax = 25/wpml;                % Maximun absorbtion

a = 0.5;

nt = 3;

hs = zeros(1,nt);
error = zeros(4,nt);
for ii = 1:nt
    ii
    tic;
    h = h/2;
    hs(ii) = h;
[node,elem] = squaremesh([-a,a,-a,a],h);
% [node,elem] = squaremesh([-.4,.6,-.4,.6],h);
Nray = 1; N = size(node,1); Ndof = N;

%% Exact ray information
xx = node(:,1)-xs;  yy = node(:,2)-ys;
rr = sqrt(xx.^2 + yy.^2);
ray = atan2(yy,xx);
ray = 0*exp(1i*ray).*(rr>10*eps);

%% Assembling
[A] = assemble_Helmholtz_matrix_with_ray_1(node,elem,omega,wpml,sigmaMax,speed,ray,fquadorder);
b = assemble_RHS_with_sing_RayFEM(node,elem,epsilon,omega,speed,xs,ys,ray,fquadorder);

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


% xx = node(:,1)-xs;  yy = node(:,2)-ys;
% rr = sqrt(xx.^2 + yy.^2);
% ub = 1i/4*besselh(0,1,omega*rr);
% aa = epsilon; bb = 2*epsilon;
% x_eps = cutoff(aa,bb,node,xs,ys);
% v = (1-x_eps).*ub;
% v(rr<aa)=0;
v = exu(node);

du = u-v;
x = node(:,1); y = node(:,2);
du(x>=max(x)-wpml)=0; du(x<= min(x)+wpml) = 0;
du(y>=max(y)-wpml)=0; du(y<= min(y)+wpml) = 0;


% idx = find( (x> min(x)+wpml).*(x<max(x)-wpml)...
%     .*(y> min(y)+wpml).*(y<max(y)-wpml) );
% 
% du = du(idx);
% v = v(idx);



max_err = norm(du,inf);
rel_max_err = norm(du,inf)/norm(v,inf);
l2_err = norm(du)*h;
rel_l2_err = norm(du)/norm(v);

% rhs = sing_rhs(epsilon,omega,node,xs,ys);
% rhsL2 = norm(rhs)*h;
error(:,ii) = [max_err,rel_max_err,l2_err,rel_l2_err]';%,rhsL2]';
toc;
end

subplot(2,2,1);
FJ_showrate(hs,error(1,:));
subplot(2,2,2);
FJ_showrate(hs,error(2,:));
subplot(2,2,3);
FJ_showrate(hs,error(3,:));
subplot(2,2,4);
FJ_showrate(hs,error(4,:));



