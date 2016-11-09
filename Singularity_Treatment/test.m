%% test for the assembling matrix and right hand side
xs = 0.1;   ys = 0.1;                  % point source location
speed = @(x) ones(size(x,1),1);    % medium speed

plt = 0;                           % plot solution or not
fquadorder = 3;                    % numerical quadrature order
pde = [];

omega = 160*pi;    
wavelength = 2*pi/omega; 
epsilon = 0.11;
NPW = 6;
h = wavelength/NPW;
h = 1/round(1/h); 
   
wpml = 2*wavelength;                 % width of PML
sigmaMax = 25/wpml;                % Maximun absorbtion

ds = 0.5;
[node,elem] = squaremesh([-ds,ds,-ds,ds],h);
% [node,elem] = squaremesh([-.6,.4,-.6,.4],h);

xx = node(:,1)-xs;  yy = node(:,2)-ys;
rr = sqrt(xx.^2 + yy.^2);
ub = 1i/4*besselh(0,1,omega*rr);
ub(rr<h/2) = 0;


figure(14);
showsolution(node,elem,real(ub));





% r1 = sort(rr);
% u1 = 1i/4*besselh(0,1,omega*r1);
% 
% xx = node(:,1);  yy = node(:,2);
% r2 = sqrt(xx.^2 + yy.^2);
% r2 = sort(r2);
% u2 = 1i/4*besselh(0,1,omega*r2);
% 
% du = u1-u2;
% norm(du(2:round(end/2)))

