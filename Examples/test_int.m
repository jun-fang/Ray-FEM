%% Test for numerical intergral
fquadorder = 9;                  % numerical quadrature order
omega = 20*pi;              % high frequency
wl = 2*pi/omega;
wpml = 2*round(wl*128)/128;%25/128;
sigmaMax = 25/wpml;


xs = 0; ys = 0;
speed = @(x) ones(size(x(:,1)));
% speed = @(x) 4/3*( 1-1/2*exp( -32*(x(:,1)-xs).^2 + (x(:,2)-ys).^2 ) );
nwl = 1/50;
sigma = nwl*wl;
source = @(x) 1/(2*pi*sigma^2)*exp(-( (x(:,1)-xs).^2 + (x(:,2)-ys).^2 )/2/sigma^2);

h = 1/512;
rh = 1/2048;
a = 1/2;

fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('Gaussian point source parameter sigma = 1/%.0f wavelength\n', wl/sigma);


%% Exact Ray-FEM

[node,elem] = squaremesh([-a,a,-a,a],h);

fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('Ray-FEM:\n omega/(2*pi) = %.2d,   1/h = %d   NPW = %.2d \n',omega/(2*pi), 1/h, wl/h);


%% Exact Ray-FEM: original right hand side
tic;
ray_ang = ex_ray(node,xs,ys);
d = sqrt((node(:,1)-xs).^2 + (node(:,2)-ys).^2);
ray = exp(1i*ray_ang).*(d>10*eps);
A = assemble_Helmholtz_matrix_with_ray_1(node,elem,omega,wpml,sigmaMax,speed,ray,fquadorder);
b = assemble_RHS_with_ray_1(node,elem,omega,source,speed,ray,fquadorder);

toc;
[~,~,isBdNode] = findboundary(elem);
freeNode = find(~isBdNode);
v = zeros(size(node(:,1)));
tic;
v(freeNode) = A(freeNode,freeNode)\b(freeNode);
toc;

% reconstruct the solution
grad = ray(:);
grad = [real(grad),imag(grad)];
temp = grad(:,1).*node(:,1) + grad(:,2).*node(:,2);

c = speed(node);    % medium speed
k = omega./c;           % wavenumber

u = v.*exp(1i*k(:).*temp);
u = sum(u,2);



%% Ray-FEM: more accurate right hand side 
fb = b;
idx = find(d<3*sigma);
for i = 1:length(idx)
    ii = idx(i);
    xi = node(ii,1); yi = node(ii,2);
    fb(ii) = RHS_integral_with_ray(xi,yi,h,rh,omega,ray(ii),source,fquadorder);
end
v(freeNode) = A(freeNode,freeNode)\fb(freeNode);
fu = v.*exp(1i*k(:).*temp);


%% S-FEM with same mesh size
sA = assemble_Helmholtz_matrix_SFEM(node,elem,omega,wpml,sigmaMax,speed,fquadorder);
sb = assemble_RHS(node,elem,source,fquadorder);
su = zeros(size(node(:,1)));
su(freeNode) = sA(freeNode,freeNode)\sb(freeNode);



%% Reference solution

fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('Reference solution S-FEM:\n omega/(2*pi) = %.2d,   1/rh = %d   NPW = %.2d \n',omega/(2*pi), 1/rh, wl/rh);

load('ex1_reference_solution_20pi.mat');

% [rnode,relem] = squaremesh([-a,a,-a,a],rh);
% tic;
% rA = assemble_Helmholtz_matrix_SFEM(rnode,relem,omega,wpml,sigmaMax,speed,fquadorder);
% rb = assemble_RHS(rnode,relem,source,fquadorder);
% toc;
% 
% [~,~,isBdNode] = findboundary(relem);
% freeNode = find(~isBdNode);
% ru = zeros(size(rnode(:,1)));
% tic;
% ru(freeNode) = rA(freeNode,freeNode)\rb(freeNode);
% toc;


%% Plot
figure(20);   % wave-field
yy = 0.25;
yy = h*ceil((yy+a)/h);

% reference solution
rn = round(sqrt(length(ru)));
ref_u = reshape(ru,rn,rn);
ryn = round(yy/rh) + 1;
rxx = -a:rh:a;
ruh = ref_u(ryn,:);

plot(rxx,real(ruh),'k');
hold on;

% solution 1
n = round(sqrt(size(node,1)));
yn = round(yy/h) + 1;
xx = -a:h:a;
uh = reshape(fu,n,n);
uh = uh(yn,:);
plot(xx,real(uh),'c^-');
hold on;
uh = reshape(u,n,n);
uh = uh(yn,:);
plot(xx,real(uh),'ro-');
hold on;
uh = reshape(su,n,n);
uh = uh(yn,:);
plot(xx,real(uh),'b+:');
hold on;

xlabel('x');
ylabel('Real part wavefield');
legend('Reference solution','Ray-FEM solution with original RHS','Ray-FEM with more accurate RHS','S-FEM solution','LOCATION','Best');
title(['Wavefield at y = ' num2str(yy-a)],'FontSize', 14)











