%% Caustics:  iterative idea

% xs = 0.5;   ys = 0.2;             % point source location
% speed = @(x) ( 1+ 0.2*sin(3*pi*(x(:,1) + 0.05)).*sin(0.5*pi*x(:,2)) );


xs = -0.2;   ys = -0.3;             % point source location
speed = @(x) ( 1+ 0.5*sin(2*pi*x(:,1)));

       

plt = 0;                           % plot solution or not
fquadorder = 4;                    % numerical quadrature order
Nray = 1;
Rest = 1;
pde = [];
pct = 1/5;
data = 'num';
opt = 1;

high_omega = 80*pi;
low_omega = sqrt(high_omega);
cmin = 1/2;
NPW = 10;

h = 1/(NPW*round(high_omega/(2*pi*cmin)));
ch = 1/(2*NPW*round(low_omega/(2*pi*cmin*2)));

wpml = 2*(2*pi*cmin)/high_omega;                 % width of PML
sigmaMax = 25/wpml;                % Maximun absorbtion
r = 4*wpml;

% a = 1;
% b = 2;
a = 1/2;
b = 1/2;
low_r = 0.5;
high_r = 0.3;


fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('omega/(2*pi) = %d,   1/h = %d   1/ch = %d,  NPW = %d \n',high_omega/(2*pi), 1/h, 1/ch, NPW);


%% Step 1: Solve the Hemholtz equation with the same source but with a relative low frequency sqrt(\omega) by Standard FEM, mesh size \omega*h = constant
fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('Step1: S-FEM, low frequency\n');

lbd = low_r + high_r;
[lnode,lelem] = squaremesh([-a-lbd,a+lbd,-b-lbd,b+lbd],h);
m = round((2*a+2*lbd)/h) + 1;
n = round((2*b+2*lbd)/h) + 1;
omega = low_omega;
[u_std] = Standard_FEM_PML_PointSource(lnode,lelem,omega,wpml,sigmaMax,xs,ys,speed,fquadorder,plt);
% FJ_showresult(lnode,lelem,real(u_std));


%% Step 2: Use NMLA to find ray directions d_c with low frequency sqrt(\omega)
fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('\nStep2: NMLA, low frequency\n');

[ux,uy] = num_derivative(u_std,h,2);

mbd = high_r;
[mnode,melem] = squaremesh([-a-mbd,a+mbd,-b-mbd,b+mbd],h);
[cnode,celem] = squaremesh([-a-mbd,a+mbd,-b-mbd,b+mbd],ch);
cN = size(cnode,1);
cnumray = zeros(cN,Nray);
cray = ex_ray_angle(cnode,xs,ys);

fprintf('NMLA time: \n');
tic;
for i = 1:cN
    x0 = cnode(i,1);  y0 = cnode(i,2);
    d0 = sqrt((x0-xs)^2 + (y0-ys)^2);
    if d0 <= r
        cnumray(i,:) = cray(i,:);
    else
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

% cdiffang = angle_error(cnumray,cray);
% norm(cdiffang,2)/norm(cray)
% norm(cdiffang,inf)

cnumray = exp(1i*cnumray);
numray1 = interpolation(cnode,celem,mnode,cnumray);

ray = ex_ray_angle(mnode,xs,ys);
ray = exp(1i*ray);
md = sqrt((mnode(:,1)-xs).^2 + (mnode(:,2)-ys).^2);
ray = ray.*(1 - (md<eps));

numray1 = numray1.*(md>r) + ray.*(md<=r);
% diffray1 = numray1 - ray;
% figure(1);
% FJ_showresult(mnode,melem,real(diffray1));

%% Step 3: Solve the original Helmholtz equation by Ray-based FEM with ray directions d_c
fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('\nStep3: Ray-FEM, high frequency \n');

omega = high_omega;
[u1] = Ray_FEM_PML_1_PointSource(mnode,melem,omega,wpml,sigmaMax,xs,ys,speed,numray1,fquadorder,plt);
FJ_showresult(mnode,melem,real(u1));

%% Step 4: NMLA to find original ray directions d_o with wavenumber k
fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('\nStep4: NMLA, high frequency \n');

[node,elem] = squaremesh([-a,a,-b,b],h);
[cnode,celem] = squaremesh([-a,a,-b,b],ch);
cN = size(cnode,1);
cnumray = zeros(cN,Nray);
cray = ex_ray_angle(cnode,xs,ys);


m = round((2*a+2*mbd)/h) + 1;
n = round((2*b+2*mbd)/h) + 1;
[ux,uy] = num_derivative(u1,h,2);

fprintf('NMLA time: \n');
tic;
for i = 1:cN
    x0 = cnode(i,1);  y0 = cnode(i,2);
    d0 = sqrt((x0-xs)^2 + (y0-ys)^2);
    if d0 <= r
        cnumray(i,:) = cray(i,:);
    else
        Rest = d0;
        if d0 <= 2*r
            Rest = 2*d0;
        end
        c0 = speed(cnode(i,:));
        [cnumray(i,:)] = NMLA_2D_2nd(x0,y0,c0,omega,Rest,mnode,melem,u1,ux,uy,pde,pct,Nray,data,opt,plt);
    end
end
toc;

% clear mnode melem;

% cdiffang = angle_error(cnumray,cray);
% norm(cdiffang,2)/norm(cray)
% norm(cdiffang,inf)

cnumray = exp(1i*cnumray);
numray2 = interpolation(cnode,celem,node,cnumray);

ray = ex_ray_angle(node,xs,ys);
ray = exp(1i*ray);
d = sqrt((node(:,1)-xs).^2 + (node(:,2)-ys).^2);
ray = ray.*(1 - (d < eps) );

numray2 = numray2.*(d>r) + ray.*(d<=r);
% diffray2 = numray2 - ray;
% figure(2);
% FJ_showresult(node,elem,real(diffray2));

% numray_dir = [real(numray2), imag(numray2)];
% numray_dir = atan2(numray_dir(:,2), numray_dir(:,1));
% numray_dir = numray_dir + 2*pi*(numray_dir<0);
% 
% diffang = angle_error(numray_dir,ex_ray_angle(node,xs,ys));
% diffang = diffang.*(d>r);
% figure(3);
% FJ_showresult(node,elem,real(diffang));
% title('NMLA angle error');


%% Step 5: Solve the original Helmholtz equation by Ray-based FEM with ray directions d_o
fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('\nStep5: Ray-FEM, high frequency \n');

omega = high_omega;
[u2] = Ray_FEM_PML_1_PointSource(node,elem,omega,wpml,sigmaMax,xs,ys,speed,numray2,fquadorder,plt);
% figure(5);
% FJ_showresult(node,elem,real(u2));
% title('Ray-FEM solution');

[X,Y] = meshgrid(-a:h:a,-b:h:b);
[m,n] = size(X);
uh = reshape(u2,m,n);

% figure(4)
% contour(X,Y,real(uh));
% title('Level set of Ray-FEM solution u_{ray}');



% ss = speed(node);
% figure(1);
% FJ_showresult(node,elem,ss);
% title('Wave speed field');

figure(10);
showsolution(node,elem,speed(node),2);
colorbar;
set(gca, 'FontSize', 18); 
set(gcf, 'Position', [250 250 600 500])

figure(20);
showsolution(node,elem,real(u2),2);
colorbar;
set(gca, 'FontSize', 18); 
set(gcf, 'Position', [1000 250 600 500])













if (0)
re = real(numray2);
im = imag(numray2);

[cX,cY] = meshgrid(-a:ch:a,-b:ch:b);
[cm,cn] = size(cX);

cre = interpolation(node,elem,cnode,re);
cim = interpolation(node,elem,cnode,im);

cre = 1/8*reshape(cre,cm,cn);
cim = 1/8*reshape(cim,cm,cn);

figure(2);
quiver(cX,cY,cre,cim)
title('Ray direction field');


figure(3);
contour(X,Y,real(uh));
title('Ray direction field');
hold on
quiver(cX,cY,cre,cim)




%% Plotting

[su] = Standard_FEM_PML_PointSource(node,elem,omega,wpml,sigmaMax,xs,ys,speed,fquadorder,plt);
N = size(node,1);
n = round(sqrt(N));
su = reshape(su,n,n);


load('test10_reference_solution_30.mat');
rn = round(1/rh) + 1;
ru = reshape(ru,rn,rn);

figure(7);
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







%% Step 6: Compute the reference solution
if (0)
NPW = 30;
rh = 1/(NPW*round(omega/(2*pi*cmin)));
[rnode,relem] = squaremesh([-a,a,-a,a],rh);
[ru] = Standard_FEM_PML_PointSource(rnode,relem,omega,wpml,sigmaMax,xs,ys,speed,fquadorder,plt);
figure(6);
FJ_showresult(rnode,relem,real(ru));
title('Reference solution');
clear rnode relem;
end
