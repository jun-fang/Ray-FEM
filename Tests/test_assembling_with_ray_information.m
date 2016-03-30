%% A test example for Ray-FEM: assembling with ray information 

omega = 100*pi;
xc = 1/10; yc = 1/10;
sigma = 1/100;

%% Load source data
source = @(x) -4*sqrt(omega)*1i*1/(2*pi*sigma^2)*exp( -( (x(:,1)-xc).^2 + (x(:,2)-yc).^2  )/(2*sigma^2) );
speed = @(x) ones(size(x(:,1)));

%% Set up
plt = 0;                   % show solution or not
fquadorder = 6;            % numerical quadrature order
NPW = 10;
h = 1/round((NPW*omega)/(2*pi));

a = 1/2;
[node,elem] = squaremesh([-a,a,-a,a],h);

wavelength = 2*pi/omega;  
wpml = 2*wavelength;               % width of PML
sigmaMax = 25/wpml;                % Maximun absorbtion

%% Exact ray information
xx = node(:,1)-xc; yy = node(:,2)-yc;
ray = atan2(yy,xx);
ray = exp(1i*ray);
ray = ray.*(1 - (abs(xx)<eps).*(abs(yy)<eps));


N = size(node,1);       % number of grid points
Nray = size(ray,2);     % number of rays crossing at each grid node
Ndof = N*Nray; 


%% Assembling the linear system

A = assemble_Helmholtz_matrix_with_ray(node,elem,omega,wpml,sigmaMax,speed,ray,fquadorder);
b = assemble_RHS_with_ray(node,elem,omega,source,speed,ray,fquadorder);

%% Boundaries
[~,~,isBdNode] = findboundary(elem);
rep_isBdNode = repmat(isBdNode,1,Nray);
isBdNode = rep_isBdNode(:);
freeNode = find(~isBdNode);


%% Solve Av=b and reconstruct the solution
v = zeros(Ndof,1);
v(freeNode) = A(freeNode,freeNode)\b(freeNode);

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


%% Show result
showresult(node,elem,real(u))


