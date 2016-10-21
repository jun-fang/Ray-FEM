%% One point source inside domain (heterogeneous medium):  iterative idea

xs = 1/10;   ys = 1/10;           % point source location
speed = @(x) (3 - 2.5*exp( -((x(:,1)+1/8).^2 + (x(:,2)-0.1).^2)/0.8^2 )) ;    % medium speed
cmin = 1/2;                      % minmum speed in the computational domain

plt = 0;                         % plot solution or not
fquadorder = 3;                  % numerical quadrature order
Nray = 1;                        % no ray crosssing, number of ray direction is 1
Rest = 1;
pde = [];
pct = 1/5;                       % percentage to specify the peak in NMLA
data = 'num';                    % use numerically computed data in NMLA 
opt = 1;                         % apply second order correction for NMLA

high_omega = 80*pi;              % high frequency
low_omega = sqrt(high_omega);    % low frequency
NPW = 10;                        % number of grid points per wavelength

h = 1/(NPW*round(high_omega/(2*pi*cmin)));    %% 1/h should be a mutiple of 40!!
ch = 1/(10*round(low_omega/(2*pi*cmin)));
  
wpml = 2*(2*pi*cmin)/high_omega;     % width of PML
sigmaMax = 25/wpml;                  % Maximun absorbtion
r = 2*wpml;                          % radius of the circle near source point

lg_a = 1;                       % large computational domain [-lg_a, lg_a]^2
md_a = 0.65;                    % middle computational domain [-md_a, md_a]^2
sm_a = 1/2;                     % physical domain [-sm_a, sm_a]^2

fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('omega/(2*pi) = %d,   1/h = %d   1/ch = %d,  NPW = %d \n',high_omega/(2*pi), 1/h, 1/ch, NPW);


%% Step 1: Solve the Hemholtz equation with the same source but with a relative low frequency sqrt(\omega) by Standard FEM, mesh size \omega*h = constant
fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('Step1: S-FEM, low frequency\n');

a = lg_a;                      % large computational domain [-lg_a, lg_a]^2
[lnode,lelem] = squaremesh([-a,a,-a,a],h);
omega = low_omega;
[u_std] = Standard_FEM_PML_PointSource(lnode,lelem,omega,wpml,sigmaMax,xs,ys,speed,fquadorder,plt);

if plt
    figure(10);
    FJ_showresult(lnode, lelem, real(u_std))
end
%% Step 2: Use NMLA to find ray directions d_c with low frequency sqrt(\omega)
fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('\nStep2: NMLA, low frequency\n');

[ux,uy] = num_derivative(u_std,h,2);    % compute derivatives numerically

a = md_a;                     % middle computational domain [-md_a, md_a]^2
[mnode,melem] = squaremesh([-a,a,-a,a],h);
[cnode,celem] = squaremesh([-a,a,-a,a],ch);
cN = size(cnode,1);
cnumray = zeros(cN,Nray);
cray = ex_ray(cnode,xs,ys,0);

fprintf('NMLA time: \n');
tic;
for i = 1:cN
    x0 = cnode(i,1);  y0 = cnode(i,2);
    d0 = sqrt((x0-xs)^2 + (y0-ys)^2);
    if d0 <= r               % ray directions near source point
        cnumray(i,:) = cray(i,:);
    else                     % ray directions computed by NMLA far away from source point
        Rest = d0;
        if d0 <= 2*r
            Rest = 2*d0;
        end
        c0 = speed(cnode(i,:));
        [cnumray(i,:)] = NMLA_2D_2nd(x0,y0,c0,omega,Rest,lnode,lelem,u_std,ux,uy,pde,pct,Nray,data,opt,plt);
    end
end
toc;

clear lnode lelem;

cnumray = exp(1i*cnumray);
numray1 = interpolation(cnode,celem,mnode,cnumray);

ray = ex_ray(mnode,xs,ys,0);
ray = exp(1i*ray);
md = sqrt((mnode(:,1)-xs).^2 + (mnode(:,2)-ys).^2);
ray = ray.*(1 - (md<eps));

% ray directions near source point
numray1 = numray1.*(md>r) + ray.*(md<=r);


%% Step 3: Solve the original Helmholtz equation by Ray-based FEM with ray directions d_c
fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('\nStep3: Ray-FEM, high frequency \n');

omega = high_omega;
[u1] = Ray_FEM_PML_1_PointSource(mnode,melem,omega,wpml,sigmaMax,xs,ys,speed,numray1,fquadorder,plt);


%% Step 4: NMLA to find original ray directions d_o with wavenumber k
fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('\nStep4: NMLA, high frequency \n');

a = sm_a;                   % physical domain [-sm_a, sm_a]^2
[node,elem] = squaremesh([-a,a,-a,a],h);
[cnode,celem] = squaremesh([-a,a,-a,a],ch);
cN = size(cnode,1);
cnumray = zeros(cN,Nray);
cray = ex_ray(cnode,xs,ys,0);

[ux,uy] = num_derivative(u1,h,2);

fprintf('NMLA time: \n');
tic;
for i = 1:cN
    i
    x0 = cnode(i,1);  y0 = cnode(i,2);
    d0 = sqrt((x0-xs)^2 + (y0-ys)^2);
    if d0 <= r              % ray directions near source point
        cnumray(i,:) = cray(i,:);
    else                    % ray directions computed by NMLA far away from source point
        Rest = d0;
        if d0 <= 2*r
            Rest = 2*d0;
        end
        c0 = speed(cnode(i,:));
        [cnumray(i,:)] = NMLA_2D_2nd(x0,y0,c0,omega,Rest,mnode,melem,u1,ux,uy,pde,pct,Nray,data,opt,plt);
    end
end
toc;

clear mnode melem;
diffcang = cnumray - cray;
cnumray = exp(1i*cnumray);
cray = exp(1i*cray);
diffcray = abs(cnumray - cray);
numray2 = interpolation(cnode,celem,node,cnumray);

ray = ex_ray(node,xs,ys,0);
ray = exp(1i*ray);
d = sqrt((node(:,1)-xs).^2 + (node(:,2)-ys).^2);
ray = ray.*(1 - (d < eps) );

% ray directions near source point
numray2 = numray2.*(d>r) + ray.*(d<=r);


%% Step 5: Solve the original Helmholtz equation by Ray-based FEM with ray directions d_o
fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('\nStep5: Ray-FEM, high frequency \n');

omega = high_omega;
[u2] = Ray_FEM_PML_1_PointSource(node,elem,omega,wpml,sigmaMax,xs,ys,speed,numray2,fquadorder,plt);

[X,Y] = meshgrid(-a:h:a,-a:h:a);
[m,n] = size(X);
uh = reshape(u2,m,n);


%% Plotting
figure(10);
showsolution(node,elem,real(u2),2);
colorbar;
set(gca, 'FontSize', 18); 
set(gcf, 'Position', [0 0 600 500])


fig2 = figure(20);
ray_field(numray2,node,20,1/5000);
axis([-0.6,0.6,-0.6,0.6]);
set(gca, 'FontSize', 18); 
set(gcf, 'Position', [0 0 550 500])


if (0)
[su] = Standard_FEM_PML_PointSource(node,elem,omega,wpml,sigmaMax,xs,ys,speed,fquadorder,plt);
N = size(node,1);
n = round(sqrt(N));
su = reshape(su,n,n);

%%
load('test9_reference_solution_40.mat');
rn = round(1/rh) + 1;
ru = reshape(ru,rn,rn);

hold off;
dy = 0.9;
yn = round(dy/h) + 1;
xx = X(yn,:);
uu = uh(yn,:);
ryn = round(dy/rh) + 1;
rxx = -0.5:rh:0.5;
ruu = ru(ryn,:);
plot(xx,real(uu),'ro-');
hold on;
ss = su(yn,:);
plot(xx,real(ss),'bs-');
hold on;
plot(rxx,real(ruu),'k');
xlabel('x');
ylabel('Wavefield');
legend('Ray-FEM solution','standard FEM solution','Reference solution','LOCATION','Best');


end




if (0)
figure(4)
contour(X,Y,real(uh));
title('Level set of Ray-FEM solution u_{ray}');

ss = speed(node);
figure(1);
FJ_showresult(node,elem,ss);
title('Wave speed field');

re = real(numray2);
im = imag(numray2);

[cX,cY] = meshgrid(-a:ch:a,-a:ch:a);
[cm,cn] = size(cX);

cre = interpolation(node,elem,cnode,re);
cim = interpolation(node,elem,cnode,im);

cre = 1/8*reshape(cre,cm,cn);
cim = 1/8*reshape(cim,cm,cn);

figure(2);
title('Ray direction field');
quiver(cX,cY,cre,cim)


figure(3);
contour(X,Y,real(uh));
title('Ray direction field');
hold on
quiver(cX,cY,cre,cim)
end




if (0)
%% Standard FEM
[su] = Standard_FEM_PML_PointSource(node,elem,omega,wpml,sigmaMax,xs,ys,speed,fquadorder,plt);
N = size(node,1);
n = round(sqrt(N));
su = reshape(su,n,n);


%% Step 6: Compute the reference solution
NPW = 30;
rh = 1/(NPW*round(omega/(2*pi*cmin)));
1/rh
[rnode,relem] = squaremesh([-a,a,-a,a],rh);
[ru] = Standard_FEM_PML_PointSource(rnode,relem,omega,wpml,sigmaMax,xs,ys,speed,fquadorder,plt);
rN = size(rnode,1);
rn = round(sqrt(rN));
ru = reshape(ru,rn,rn);
% save('reference_solutions9.mat','omega','rh','ru');

%%
hold off;
dy = 0.9;
yn = round(dy/h) + 1;
xx = X(yn,:);
uu = uh(yn,:);
suu = su(yn,:);
ryn = round(dy/rh) + 1;
rxx = -0.5:rh:0.5;
ruu = ru(ryn,:);
plot(xx,real(uu),'ro-');
hold on;
plot(xx,real(suu),'b^-');
hold on;
plot(rxx,real(ruu),'k');
xlabel('x');
ylabel('Wavefield');
legend('Ray-FEM solution','Standard FEM solution','Reference solution','LOCATION','Best');
% axis equal;
end