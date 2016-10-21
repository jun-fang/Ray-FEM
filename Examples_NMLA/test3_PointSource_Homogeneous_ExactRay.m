%% One point source inside domain: homogeneous medium

% point source information
xs = -0.4;   ys = -0.4;             % point source location
speed = @(x) ones(size(x,1),1);     % medium speed

% parameters for numerical methods
plt = 0;                            % plot solution or not
fquadorder = 4;                     % numerical quadrature order

% physical parameters
high_omega = 80*pi;                % high frequency
NPW = 6;                            % number of grid points per wavelength
omega = high_omega;

h = 1/(20*round(high_omega*NPW/(2*pi*20)));             % fine mesh

wl = 2*pi/high_omega;               % wavelength 
wpml = h*ceil(wl/h);                % width of PML
sigmaMax = 25/wpml;                 % Maximun absorbtion
r = 8*wpml;

sigma = h/2;
source = @(x) -4*sqrt(omega)*1i*1/(2*pi*sigma^2)*exp( -( (x(:,1)-xs).^2 + (x(:,2)-ys).^2  )/(2*sigma^2) );
a = 1/2;
[node,elem] = squaremesh([-a,a,-a,a],h);
N = size(node,1);
Ndof = N; Nray = 1;

%% Exact ray information
xx = node(:,1)-xs; yy = node(:,2)-ys;
ray = atan2(yy,xx);
ray = exp(1i*ray);
ray = ray.*(1 - (abs(xx)<eps).*(abs(yy)<eps));

A = assemble_Helmholtz_matrix_with_ray_1(node,elem,omega,wpml,sigmaMax,speed,ray,fquadorder);
b = assemble_RHS_with_ray_1(node,elem,omega,source,speed,ray,fquadorder);
b = b/(h*h);

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


% %% Show result
% figure(10);
% hold off;
% subplot(2,1,1);
% nr = 30;
% r1 = wl*(nr+0.6);
% r2 = wl*(nr+2.6);
% theta1 = 0;
% theta2 = pi/2;
% theta = theta1:pi/100:theta2;
% r = r1:1/1000:r2;
% x1 = r1.*cos(theta) + xs;
% y1 = r1.*sin(theta) + ys;
% x2 = r2.*cos(theta) + xs;
% y2 = r2.*sin(theta) + ys;
% x3 = r.*cos(theta1) + xs;
% y3 = r.*sin(theta1) + ys;
% x4 = r.*cos(theta2) + xs;
% y4 = r.*sin(theta2) + ys;
% % plot(x1,y1,x2,y2,x3,y3,x4,y4,'wo');
% plot(x1,y1,x2,y2,x3,y3,x4,y4,'r','LineWidth',4);
% hold on;
% showsolution(node,elem,real(u),2)
% axis equal; axis tight; 
% colorbar;
% 
% subplot(2,1,2);
% polar_plot_solution(omega,node,xs,ys,r1,r2,theta1,theta2,u,'spline',0);
% axis equal; axis tight;

%%
figure(10);
hold off;
% sp1 = subplot(2,1,1);
positionVector1 = [0.15, 0.25, 0.7, 0.7];
subplot('Position',positionVector1)

showsolution(node,elem,real(u),2)
axis equal; axis off; 
colorbar;
hold on

nr = 30;
r1 = wl*(nr+0.6);
r2 = wl*(nr+2.6);
theta1 = 0;
theta2 = pi/2;
theta = theta1:pi/100:theta2;
r = r1:1/1000:r2;
x1 = r1.*cos(theta) + xs;
y1 = r1.*sin(theta) + ys;
x2 = r2.*cos(theta) + xs;
y2 = r2.*sin(theta) + ys;
x3 = r.*cos(theta1) + xs;
y3 = r.*sin(theta1) + ys;
x4 = r.*cos(theta2) + xs;
y4 = r.*sin(theta2) + ys;
p = plot(x1,y1,x2,y2,x3,y3,x4,y4,'Color','r','LineWidth',4);
axis equal; 
axis([-0.5,0.5,-0.5,0.5]);
% set(p,'Color','red');
hold all;

positionVector2 = [0.15, 0.1, 0.7, 0.1];
subplot('Position',positionVector2)

% sp2 = subplot(2,1,2);
polar_plot_solution(omega,node,xs,ys,r1,r2,theta1,theta2,u,'spline',0);
axis equal; axis tight;
% colorbar;
set(gca, 'FontSize', 18); 
% set(gcf, 'Position', [1000 250 600 500])
