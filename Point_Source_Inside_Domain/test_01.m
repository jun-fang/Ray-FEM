global omega;
omega = 100*pi;
xc = 2; yc = 2;
pde = Helmholtz_data1;
sigma = 1/100;
   
%% Load source data
% source = @(x) sqrt(omega)*besselh(0,1,omega*sqrt((x(:,1)-xc).^2 + (x(:,2)-yc).^2));
source = @(x) -4*sqrt(omega)*1i*1/(2*pi*sigma^2)*exp( -( (x(:,1)-xc).^2 + (x(:,2)-yc).^2  )/(2*sigma^2) );
speed = @(x) ones(size(x(:,1)));

%% Set up
plt = 0;                   % show solution or not
fquadorder = 4;            % numerical quadrature order
NPW = 8;
h = 1/round((NPW*omega)/(2*pi));

wavelength = 2*pi/omega;
wpml = 4*round(wavelength/h)*h;               % width of PML
sigmaMax = 50/wpml;                % Maximun absorbtion

low_omega = sqrt(omega);
bd = 1/4/sqrt(low_omega);
bd = round(bd/h)*h;

a = 1/2+wpml;
[node,elem] = squaremesh([-a,a,-a,a],h);
% [node,elem] = squaremesh_annulus([-a,a,-a,a],[-bd,bd,-bd,bd],h);


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

[u,A2,v,b2] = Ray_FEM_PML_1(node,elem,omega,wpml,sigmaMax,pde,ray,fquadorder,plt);
da = A-A2;
norm(da(:))
norm(b -b2)



