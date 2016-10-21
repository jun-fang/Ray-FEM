%% One point source inside domain (homogeneous medium)
clear;
xc = 0;   yc = 0;           % point source location
speed = @(x) ones(size(x(:,1)));
   
plt = 0;                         % plot solution or not
fquadorder = 9;                  % numerical quadrature order
Nray = 1;                        % no ray crosssing, number of ray direction is 1

omega = 20*pi;              % high frequency
NPW = 10;                         % number of grid points per wavelength
pow = round(log2((NPW*omega)/(2*pi)));
   
wl = 2*pi/omega;
npml = NPW;
h = 1/(2^pow - 2*npml);
wpml = h*npml;
sigmaMax = 25/wpml;

a = 1/2 + wpml;
% [node,elem] = squaremesh([-a,a,-a,a],h);
node = [-a,-a; a,-a; a,a; -a,a];
elem = [2,3,1; 4,1,3];      
for k = 1:pow
    [node,elem] = uniformbisect(node,elem);
end
h - 2*a/2^pow    


fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('Ray-FEM:\n omega/(2*pi) = %.2d,   1/h = %d   NPW = %.2d \n',omega/(2*pi), 1/h, wl/h);

% sigma = 1/100;
% source = @(x) -4*sqrt(omega)*1i*1/(2*pi*sigma^2)*exp( -( (x(:,1)-xc).^2 + (x(:,2)-yc).^2  )/(2*sigma^2) );


%% Exact Ray-FEM
ray_ang = ex_ray(node,xc,yc);
d = sqrt((node(:,1)-xc).^2 + (node(:,2)-yc).^2);
ray = exp(1i*ray_ang).*(d>10*eps);
[A,M] = assemble_Helmholtz_matrix_with_ray_1(node,elem,omega,wpml,sigmaMax,speed,ray,fquadorder);

% [A,M,D] = assemble_Helmholtz_matrix_SFEM(node,elem,omega,wpml,sigmaMax,speed,fquadorder);

% n = round(2*a/h + 1);
% xn = round((xc + a)/h);
% yn = round((yc + a)/h);
% xyn = n*xn + yn + 1;  

xyn = find(d<10*eps);
b = M(:,xyn)/(1/2*h^2);

% b = assemble_ref_RHS(node,elem,h,xc,yc,9);
% b = b/(1/2*h*h);
% norm(b1-b)
% b = assemble_RHS_with_ray_1(node,elem,omega,source,speed,ray,fquadorder);


[~,~,isBdNode] = findboundary(elem);
freeNode = find(~isBdNode);
N = size(node,1);
v = zeros(N,1);
tic;
v(freeNode) = A(freeNode,freeNode)\b(freeNode);
toc;


%% reconstruct the solution
grad = ray(:);
grad = [real(grad),imag(grad)];
repnode = repmat(node,Nray,1);
temp = grad(:,1).*repnode(:,1) + grad(:,2).*repnode(:,2);

c = speed(node);    % medium speed
k = omega./c;           % wavenumber
kk = repmat(k,1,Nray);

u = v.*exp(1i*kk(:).*temp);
u = reshape(u,N,Nray);
u = sum(u,2);
% figure(1);
% FJ_showresult(node,elem,real(u));


%% S-FEM with same mesh
sA = assemble_Helmholtz_matrix_SFEM(node,elem,omega,wpml,sigmaMax,speed,fquadorder);
sb = assemble_ref_RHS(node,elem,h,xc,yc,9);
sb = sb/(1/2*h*h);
sv = v;
sv(freeNode) = sA(freeNode,freeNode)\sb(freeNode);



%% reference solution? S-FEM
rpow = 4;
rnode = node; relem = elem;
for i = 1:rpow
    [rnode,relem] = uniformbisect(rnode,relem);
end
rh = h/2^rpow; 
% rh = 1/1600;

fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('Reference S-FEM:\n omega/(2*pi) = %.2d,   1/rh = %d   NPW = %.2d \n',omega/(2*pi), 1/rh, wl/rh);
% [rnode,relem] = squaremesh([-a,a,-a,a],rh);
% [rnode,relem] = uniformbisect(rnode,relem);
rA = assemble_Helmholtz_matrix_SFEM(rnode,relem,omega,wpml,sigmaMax,speed,fquadorder);

% rn = round(2*a/rh + 1);
% rxn = round((xc + a)/rh);
% ryn = round((yc + a)/rh);
% rxyn = rn*rxn + ryn + 1; 
% rb1 = rM(:,rxyn)/(1/2*rh^2);

rb = assemble_ref_RHS(rnode,relem,h,xc,yc,9);
rb = rb/(1/2*h*h);
% rb = assemble_RHS(rnode,relem,source,fquadorder);

[~,~,isBdNode] = findboundary(relem);
freeNode = find(~isBdNode);
ru = zeros(size(rnode(:,1)));
tic;
ru(freeNode) = rA(freeNode,freeNode)\rb(freeNode);
toc;
% figure(2);
% FJ_showresult(rnode,relem,real(ru));


%% Plotting
% figure(1);
% FJ_showresult(node,elem,real(u));
% 
% figure(2);
% FJ_showresult(rnode,relem,real(ru));
[nnode,u] = transform_node(node,u);
[nnnode,rru] = transform_node(rnode,ru);
[snode,su] = transform_node(node,sv);

figure(2)
show_ray_solution(a,h,rh,u,su,rru,0.3)


%% Test for approximate delta function




