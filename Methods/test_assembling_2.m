%% A test example for Ray-FEM: assembling with ray information 

omega = 40*pi;
xc = 1/10; yc = 1/10;
sigma = 1/100;

%% Load source data
source = @(x) -4*sqrt(omega)*1i*1/(2*pi*sigma^2)*exp( -( (x(:,1)-xc).^2 + (x(:,2)-yc).^2  )/(2*sigma^2) );
speed = @(x) ones(size(x(:,1)));

%% Set up
plt = 0;                   % show solution or not
fquadorder = 3;            % numerical quadrature order
NPW = 8;
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
exray = ray.*(1 - (abs(xx)<eps).*(abs(yy)<eps));


N = size(node,1);       % number of grid points
Nray = size(ray,2);     % number of rays crossing at each grid node

ray = cell(N,1);
for i = 1:N
    ray{i} = exray(i);
end

ray_num = zeros(N,1);     % number of rays at each grid point
ray_dof = zeros(N,1);     % ray_dof(n) = sum(ray_num(1:n)) 

temp = 0;
for n = 1:N
    ray_num(n) = size(ray{n},2);
    ray_dof(n) = temp + ray_num(n);
    temp = ray_dof(n);
end

Ndof = temp;              % total degree of freedom


%% Assembling the linear system
A = assemble_Helmholtz_matrix_with_ray_2(node,elem,omega,wpml,sigmaMax,speed,ray,fquadorder);
b = assemble_RHS_with_ray_2(node,elem,omega,source,speed,ray,fquadorder);


%% Boundary conditions
[bdnode,~,~] = findboundary(elem);
isFreeNode = ones(Ndof,1);
Nbd = length(bdnode);
for nn = 1:Nbd
    nbd = bdnode(nn);
    ind = ray_dof(nbd) - ray_num(nbd) + 1:1:ray_dof(nbd);
    isFreeNode(ind) = 0;
end
freeNode = find(isFreeNode);


%% Solve the linear system Au = b
v = zeros(Ndof,1);
v(freeNode) = A(freeNode,freeNode)\b(freeNode);


%% re-construct the solution at grid nodes
c = speed(node);      % medium speed 
k = omega./c;             % wavenumber
u = zeros(N,1);
for n = 1:N
    grad = ray{n};
    grad = transpose(grad);
    grad = [real(grad), imag(grad)];    
    tempn = node(n,1)*grad(:,1) + node(n,2)*grad(:,2);
    ni = ray_dof(n) - ray_num(n);
    nii = 1:ray_num(n);
    nii = nii' + ni;
    u(n) = sum(v(nii).*exp(1i*k(n)*tempn));
end


%% Show result
showresult(node,elem,real(u))


