%% One point source inside domain (homogeneous medium)
clear;
xc = -3/8;   yc = -3/8;           % point source location
speed = @(x) ones(size(x(:,1)));
   
plt = 0;                         % plot solution or not
fquadorder = 3;                  % numerical quadrature order
Nray = 1;                        % no ray crosssing, number of ray direction is 1
    
omega = 20*pi;              % high frequency
NPW = 16;                         % number of grid points per wavelength
h = 1/2^(round(log2((NPW*omega)/(2*pi))));    

wl = 2*pi/omega;
wpml = 2*ceil(wl/h)*h;
sigmaMax = 25/wpml;

a = 1/2 + h*ceil(wpml/h);
[node,elem] = squaremesh([-a,a,-a,a],h);

fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('Ray-FEM:\n omega/(2*pi) = %.2d,   1/h = %d   NPW = %d \n',omega/(2*pi), 1/h, NPW);


sigma = 1/100;
source = @(x) -4*sqrt(omega)*1i*1/(2*pi*sigma^2)*exp( -( (x(:,1)-xc).^2 + (x(:,2)-yc).^2  )/(2*sigma^2) );


%% Ray-FEM
ray = ex_ray(node,xc,yc,1);
[A,M] = assemble_Helmholtz_matrix_with_ray_1(node,elem,omega,wpml,sigmaMax,speed,ray,fquadorder);
b = assemble_RHS_with_ray_1(node,elem,omega,source,speed,ray,fquadorder);

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


%% reference solution
NPW = 80;
rh = 1/2^(round(log2((NPW*omega)/(2*pi)))); 
% rh = 1/1600;

fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('Reference S-FEM:\n omega/(2*pi) = %.2d,   1/rh = %d   NPW = %d \n',omega/(2*pi), 1/rh, NPW);
[rnode,relem] = squaremesh([-a,a,-a,a],rh);
[rA,rM] = assemble_Helmholtz_matrix_SFEM(rnode,relem,omega,wpml,sigmaMax,speed,fquadorder);
rb = assemble_RHS(rnode,relem,source,fquadorder);
% rb = rb/(sum(b)*rh^2);
[~,~,isBdNode] = findboundary(relem);
freeNode = find(~isBdNode);
ru = zeros(size(rnode(:,1)));
ru(freeNode) = rA(freeNode,freeNode)\rb(freeNode);
% figure(2);
% FJ_showresult(rnode,relem,real(ru));


%% Plotting
show_ray_solution(a,h,rh,u,ru,0.3)





