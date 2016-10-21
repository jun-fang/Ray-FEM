%% Ray information on the computational domain or physical domain ??
% verifying the necessity of ray information on the whole computational
% domain, not just on the physical domain

xc = 1/8;   yc = 1/16;              % point source location
speed = @(x) ones(size(x(:,1)));    % medium speed
   
plt = 0;                            % plot solution or not
fquadorder = 3;                     % numerical quadrature order
Nray = 1;                           % no ray crosssing, number of ray direction is 1

omega = 20*pi;                     % frequency
sigma = 1/100;
source = @(x) -4*sqrt(omega)*1i*1/(2*pi*sigma^2)*exp( -( (x(:,1)-xc).^2 + (x(:,2)-yc).^2  )/(2*sigma^2) );
NPW = 8;                            % number of grid points per wavelength
h = 1/2^(round(log2((NPW*omega)/(2*pi))));    % mesh size

wl = 2*pi/omega;                    % wavelength
wpml = 2*h*ceil(wl/h);              % width of PML
sigmaMax = 25/wpml;                 % maxmium absorption
a = 1/2 + wpml;                     % computational domain [-a,a]^2

[node,elem] = squaremesh([-a,a,-a,a],h);
N = size(node,1);

%% ray information
xx = node(:,1) - xc; yy = node(:,2) - yc;
rayang = atan2(yy,xx);
ray = exp(1i*rayang);
ray = ray.*(norm([xx,yy])>eps);


%% solve Av=b with exact ray on computational domain [-a,a]^2
A = assemble_Helmholtz_matrix_with_ray_1(node,elem,omega,wpml,sigmaMax,speed,ray,fquadorder);
b = assemble_RHS_with_ray_1(node,elem,omega,source,speed,ray,fquadorder);

[~,~,isBdNode] = findboundary(elem);
freeNode = find(~isBdNode);
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
u1 = u;

%% solve Av=b with exact ray on physical domain [-1/2,1/2]^2
phy_dom = (node(:,1)<=1/2).*(node(:,1)>=-1/2).*(node(:,2)<=1/2).*(node(:,2)>=-1/2);
ray = ray.*phy_dom;

%% solve Av=b
A = assemble_Helmholtz_matrix_with_ray_1(node,elem,omega,wpml,sigmaMax,speed,ray,fquadorder);
b = assemble_RHS_with_ray_1(node,elem,omega,source,speed,ray,fquadorder);

[~,bdEdge,isBdNode] = findboundary(elem);
freeNode = find(~isBdNode);
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
u2 = u;


%% compare two solutions
figure(1);
FJ_showresult(node,elem,real(u1));
figure(2);
FJ_showresult(node,elem,real(u2));

du = u1-u2;
du = du(phy_dom>0);
norm(du,inf)



% %% Plotting
% [su] = Standard_FEM_PML_PointSource(node,elem,omega,wpml,sigmaMax,xs,ys,speed,fquadorder,plt);
% N = size(node,1);
% n = round(sqrt(N));
% su = reshape(su,n,n);
% 
% 
% load('test9_reference_solution_40.mat');
% rn = round(1/rh) + 1;
% ru = reshape(ru,rn,rn);
% 
% hold off;
% dy = 0.9;
% yn = round(dy/h) + 1;
% xx = X(yn,:);
% uu = uh(yn,:);
% ryn = round(dy/rh) + 1;
% rxx = -0.5:rh:0.5;
% ruu = ru(ryn,:);
% plot(xx,real(uu),'ro-');
% hold on;
% ss = su(yn,:);
% plot(xx,real(ss),'bs-');
% hold on;
% plot(rxx,real(ruu),'k');
% xlabel('x');
% ylabel('Wavefield');
% legend('Ray-FEM solution','standard FEM solution','Reference solution','LOCATION','Best');
% 
% 
% 
% 
% 
% 
% 
% if (0)
% figure(4)
% contour(X,Y,real(uh));
% title('Level set of Ray-FEM solution u_{ray}');
% 
% ss = speed(node);
% figure(1);
% FJ_showresult(node,elem,ss);
% title('Wave speed field');
% 
% re = real(numray2);
% im = imag(numray2);
% 
% [cX,cY] = meshgrid(-a:ch:a,-a:ch:a);
% [cm,cn] = size(cX);
% 
% cre = interpolation(node,elem,cnode,re);
% cim = interpolation(node,elem,cnode,im);
% 
% cre = 1/8*reshape(cre,cm,cn);
% cim = 1/8*reshape(cim,cm,cn);
% 
% figure(2);
% title('Ray direction field');
% quiver(cX,cY,cre,cim)
% 
% 
% figure(3);
% contour(X,Y,real(uh));
% title('Ray direction field');
% hold on
% quiver(cX,cY,cre,cim)
% end
% 
% 
% 
% 
% if (0)
% %% Standard FEM
% [su] = Standard_FEM_PML_PointSource(node,elem,omega,wpml,sigmaMax,xs,ys,speed,fquadorder,plt);
% N = size(node,1);
% n = round(sqrt(N));
% su = reshape(su,n,n);
% 
% 
% %% Step 6: Compute the reference solution
% NPW = 30;
% rh = 1/(NPW*round(omega/(2*pi*cmin)));
% 1/rh
% [rnode,relem] = squaremesh([-a,a,-a,a],rh);
% [ru] = Standard_FEM_PML_PointSource(rnode,relem,omega,wpml,sigmaMax,xs,ys,speed,fquadorder,plt);
% rN = size(rnode,1);
% rn = round(sqrt(rN));
% ru = reshape(ru,rn,rn);
% % save('reference_solutions9.mat','omega','rh','ru');
% 
% %%
% hold off;
% dy = 0.9;
% yn = round(dy/h) + 1;
% xx = X(yn,:);
% uu = uh(yn,:);
% suu = su(yn,:);
% ryn = round(dy/rh) + 1;
% rxx = -0.5:rh:0.5;
% ruu = ru(ryn,:);
% plot(xx,real(uu),'ro-');
% hold on;
% plot(xx,real(suu),'b^-');
% hold on;
% plot(rxx,real(ruu),'k');
% xlabel('x');
% ylabel('Wavefield');
% legend('Ray-FEM solution','Standard FEM solution','Reference solution','LOCATION','Best');
% % axis equal;
% end