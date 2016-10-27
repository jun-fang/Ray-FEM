%% Homogeneous medium
% check: 1. assembling matirx 2. right hand side 3. smoothness of function
% 4.suitable parameters


xs = 0;   ys = 0;                  % point source location
speed = @(x) ones(size(x,1),1);    % medium speed

fquadorder = 3;                    % numerical quadrature order

omega = 160*pi;    
wavelength = 2*pi/omega; 
epsilon = 0.137;%0.11/wavelength;
NPW = 6;
h = wavelength/NPW;
h = 1/round(1/h);
   
wpml = 0.08; %h*round(wavelength/h);                 % width of PML
sigmaMax = 25/wpml;                % Maximun absorbtion

a = 1/2;
[node,elem] = squaremesh([-a,a,-a,a],h);
Nray = 1; N = size(node,1); Ndof = N; n = round(sqrt(N));

xx = node(:,1)-xs;  yy = node(:,2)-ys;
rr = sqrt(xx.^2 + yy.^2);
ray = atan2(yy,xx);
ray = exp(1i*ray).*(rr>eps);


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
aa = epsilon;  bb = 2*epsilon;
x_eps = cutoff(aa,bb,node,xs,ys);
v = ub.*x_eps;
us = u;

uu = us + v;
toc;

figure(11);
showsolution(node,elem,real(us));

figure(12);
showsolution(node,elem,real(uu));

figure(13);
showsolution(node,elem,real(x_eps));

du = ub -uu;
x = node(:,1); y = node(:,2);
du(x>=a-wpml)=0; du(x<= -a+wpml) = 0;
du(y>=a-wpml)=0; du(y<= -a+wpml) = 0;
figure(14);
showsolution(node,elem,real(du));


%% plot wavefield 
figure(22);
px = 0.4;
px = round(px/h)*h;
indx = 1 + round((px+a)/h);
n = round(sqrt(N));
uh = reshape(uu,n,n);
uex = reshape(ub,n,n);

xx = -a:h:a;
yh = uh(indx,:);
yex = uex(indx,:);

plot(xx,real(yh),'ro-');
hold on;
plot(xx,real(yex));
title('real part of the wavefield at y = 0.4','FontSize', 30);
xlabel('x','FontSize', 30)
ylabel('wavefield','FontSize', 30)
legend('Numerical solution','Exact solution','Location','best');







