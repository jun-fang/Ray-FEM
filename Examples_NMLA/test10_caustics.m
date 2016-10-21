clear;   
xs = -0.2;   ys = -0.3;             % point source location
speed = @(x) ( 1+ 0.5*sin(2*pi*x(:,1)));

plt = 0;                           % plot solution or not
fquadorder = 3;                    % numerical quadrature order
Nray = 0;
Rest = 1;    
pde = [];
pct = 1/2; 
data = 'num';
opt = 0;

high_omega = 60*pi;
low_omega = 2*sqrt(high_omega);
cmin = 1/2;
NPW = 10;

h = 1/(NPW*round(high_omega/(2*pi*cmin)));
ch = 1/(2*NPW*round(low_omega/(2*pi*cmin*2)));

minwavelength = (2*pi*cmin)/high_omega;

wpml = 10*minwavelength;                 % width of PML
sigmaMax = 50/wpml;                % Maximun absorbtion
r = 4*minwavelength;

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

lbd =  low_r + high_r;
[lnode,lelem] = squaremesh([-a-lbd,a+lbd,-b-lbd,b+lbd],h);
m = round((2*a+2*lbd)/h) + 1;
n = round((2*b+2*lbd)/h) + 1;
omega = low_omega;
[u_std] = Standard_FEM_PML_PointSource(lnode,lelem,omega,wpml,sigmaMax,xs,ys,speed,fquadorder,plt);
% FJ_showresult(lnode,lelem,real(u_std));


%% Step 2: Use NMLA to find ray directions d_c with low frequency sqrt(\omega)
fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('\nStep2: NMLA, low frequency\n');

[ux,uy] = num_derivative_2(u_std,h,m,n,2);

mbd = high_r;
[mnode,melem] = squaremesh([-a-mbd,a+mbd,-b-mbd,b+mbd],h);
m = round((2*a+2*mbd)/h) + 1;
n = round((2*b+2*mbd)/h) + 1;
mN = size(mnode,1);
mnumray = cell(mN,1);
mray = ex_ray(mnode,xs,ys,0);

fprintf('NMLA time: \n');
tic;
for i = 1:mN
    i
    mN
    x0 = mnode(i,1);  y0 = mnode(i,2);
    d0 = sqrt((x0-xs)^2 + (y0-ys)^2);
    if d0 <= r
        if d0 < eps
            mnumray{i} = 0;
        else
            mnumray{i} = exp(1i*mray(i,:));
        end
    else
        Rest = d0;
        if d0 <= 2*r
            Rest = 2*d0;
        end
        c0 = speed(mnode(i,:));
        [ang] = NMLA_2D_2nd(x0,y0,c0,omega,Rest,lnode,lelem,u_std,ux,uy,pde,pct,Nray,data,opt,plt);
        mnumray{i} = exp(1i*ang);
    end
end
toc;
figure(1);
ray_field(mnumray,mnode,10);
clear lnode lelem;


%% Step 3: Solve the original Helmholtz equation by Ray-based FEM with ray directions d_c
fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('\nStep3: Ray-FEM, high frequency \n');

omega = high_omega;
[u1] = Ray_FEM_PML_22_PointSource(mnode,melem,omega,wpml,sigmaMax,xs,ys,speed,mnumray,fquadorder,plt);
figure(2);
FJ_showresult(mnode,melem,real(u1));

%% Step 4: NMLA to find original ray directions d_o with wavenumber k
fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('\nStep4: NMLA, high frequency \n');

[ux,uy] = num_derivative_2(u1,h,m,n,2);

[node,elem] = squaremesh([-a,a,-b,b],h);
N = size(node,1);
numray = cell(N,1);
ray = ex_ray(node,xs,ys,0);

fprintf('NMLA time: \n');
tic;
for i = 1:N
    i 
    N
    x0 = node(i,1);  y0 = node(i,2);
    d0 = sqrt((x0-xs)^2 + (y0-ys)^2);
    if d0 <= r
        if d0 < eps
            numray{i} = 0;
        else
            numray{i} = exp(1i*ray(i,:));
        end
    else
        Rest = d0;
        if d0 <= 2*r
            Rest = 2*d0;
        end
        c0 = speed(node(i,:));
        [ang] = NMLA_2D_2nd(x0,y0,c0,omega,Rest,mnode,melem,u1,ux,uy,pde,pct,Nray,data,opt,plt);
        numray{i} = exp(1i*ang);
    end
end
toc;
figure(3);
ray_field(numray,node,10);
clear mnode melem;


%% Step 5: Solve the original Helmholtz equation by Ray-based FEM with ray directions d_o
fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('\nStep5: Ray-FEM, high frequency \n');

omega = high_omega;
[u2] = Ray_FEM_PML_22_PointSource(node,elem,omega,wpml,sigmaMax,xs,ys,speed,numray,fquadorder,plt);
figure(4);
FJ_showresult(node,elem,real(u2));
title('Ray-FEM solution');

[X,Y] = meshgrid(-a:h:a,-b:h:b);
[m,n] = size(X);
uh = reshape(u2,m,n);

figure(5)
contour(X,Y,real(uh));
title('Level set of Ray-FEM solution u_{ray}');


%% Plotting

[su] = Standard_FEM_PML_PointSource(node,elem,omega,wpml,sigmaMax,xs,ys,speed,fquadorder,plt);
N = size(node,1);
n = round(sqrt(N));
su = reshape(su,n,n);

%% Reference
fh = 1/2000;
[fnode,felem] = squaremesh([-0.5,0.5,-0.5,0.5],fh);
[uf] = Standard_FEM_PML_PointSource(fnode,felem,omega,wpml,sigmaMax,xs,ys,speed,fquadorder,plt);

%% plot
a = 0.5;
[X,Y] = meshgrid(-a:h:a,-b:h:b);
[m,n] = size(X);

u2 = reshape(u2,m,n);
us = reshape(su,m,n);


fn = round(1/fh) + 1;
uf = reshape(uf,fn,fn);


figure(7);
hold off;
dy = 0.8;
yn = round(dy/h) + 1;
xx = X(yn,:);
uu = u2(yn,:);
ryn = round(dy/fh) + 1;
rxx = -0.5:fh:0.5;
ruu = uf(ryn,:);
plot(xx,real(uu),'ro-');
hold on;
ss = us(yn,:);
plot(xx,real(ss),'bs-');
hold on;
plot(rxx,real(ruu),'k');
xlabel('x');
ylabel('Wavefield');
legend('Ray-FEM solution','standard FEM solution','Reference solution','LOCATION','Best');
