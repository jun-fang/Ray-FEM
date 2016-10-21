%% test for the assembling matrix and right hand side
xs = 0;   ys = 0;                  % point source location
speed = @(x) ones(size(x,1),1);    % medium speed

plt = 0;                           % plot solution or not
fquadorder = 3;                    % numerical quadrature order
pde = [];

omega = 40*pi;    
wavelength = 2*pi/omega; 
epsilon = 0.11;
NPW = 100;
h = wavelength/NPW;
h = 1/round(1/h); 
   
wpml = 2*wavelength;                 % width of PML
sigmaMax = 25/wpml;                % Maximun absorbtion

ds = 1/2;
[node,elem] = squaremesh([-ds,ds,-ds,ds],h);

xx = node(:,1)-xs;  yy = node(:,2)-ys;
rr = sqrt(xx.^2 + yy.^2);
ray = atan2(yy,xx);
ray = exp(1i*ray).*(rr>eps);


% [A1] = assemble_Helmholtz_matrix_with_ray_1(node,elem,omega,wpml,sigmaMax,speed,ray,fquadorder);
% b1 = assemble_RHS_with_ray_sing(node,elem,epsilon,omega,speed,ray,fquadorder);
p = node;
a = epsilon; b = 2*a;
% cf = cutoff(a,b,p);
cl = cutoff_laplacian(a,b,node);
% cg = cutoff_gradient(a,b,node);

% figure(11);
% showsolution(node,elem,real(cf));

figure(12);
showsolution(node,elem,real(cl));

% figure(13);
% showsolution(node,elem,real(cg(:,1)));
% 
% figure(14);
% showsolution(node,elem,real(cg(:,2)));
