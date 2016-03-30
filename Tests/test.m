pde = Helmholtz_data3;
    
global omega;
global a;
a = 1/2;
omega = 100*pi;
NPW = 8;

h = 1/round((NPW*omega)/(2*pi));
wavelength = 2*pi/omega;
wl = wavelength;
T = h/8;

r = -wl: T : wl;
r = r';
est_ang = oangs1(1);
x = r*cos(est_ang);
y = r*sin(est_ang);
node = [x, y];
u = pde.ex_u(node);

% tol = 1e-3;
% L = 10;
[z] = matrixpencil(u, x)
acos(asin(max(imag(z)))/(omega*T))
% ang = [0;pi/2;pi;3*pi/2];
% exp(1i*omega*h*cos(ang))

% t = 0:24 ;
% y = sin(3*t+pi/4)+randn(size(t))/10 ;
% z = mpencil(y,2,8) ;
% log(z)








%% test for optimization
if (0)
global omega;

for i = 1: 1
    omega = 100*pi*i;
    test1_optimization(omega);
end

end


%% NMLA
if (0)
%% NMLA 
xs = -0.2; ys = -0.3;
Nray = 0;
Rest = 1;
pde = [];
pct = 1/2;      
data = 'num';
opt = 0;

minwavelength = (2*pi*cmin)/omega;
r = 4*minwavelength;

m = size(rnode,1);
m = round(sqrt(m));
n = m;

[ux,uy] = num_derivative_2(ru,rh,m,n,2);
a = 0.5;
b = 0.5;
mbd = 0;
% mbd = high_r;
% [mnode,melem] = squaremesh([-a-mbd,a+mbd,-b-mbd,b+mbd],h);
ch = 1/60;
[cnode,celem] = squaremesh([-a-mbd,a+mbd,-b-mbd,b+mbd],ch);

NPW = 6;
h = 1/(NPW*round(omega/(2*pi*cmin)));
[node,elem] = squaremesh([-0.5,0.5,-0.5,0.5],h);
N = size(node,1);
numray = cell(N,1);

cN = size(cnode,1);
cnumray = zeros(cN,1);
cray = ex_ray(cnode,xs,ys,1);
ray = ex_ray(node,xs,ys,1);


Nray = 1;
opt = 1;
fprintf('NMLA time: \n');
tic;
for i = 1:cN
    x0 = cnode(i,1);  y0 = cnode(i,2);
    d0 = sqrt((x0-xs)^2 + (y0-ys)^2);
    Rest = d0;
    c0 = speed(cnode(i,:));
    if x0 >= 0 || y0 <= 0.2
        [ang] = NMLA_2D_2nd(x0,y0,c0,omega,Rest,rnode,relem,ru,ux,uy,pde,pct,Nray,data,opt,plt);
    end
    cnumray(i) = exp(1i*ang);
end
toc;
numray1 = interpolation(cnode,celem,node,cnumray);
numray1 = ray_convert(numray1,2);

opt = 0;
tic;
for i = 1:N
    x0 = node(i,1);  y0 = node(i,2);
    d0 = sqrt((x0-xs)^2 + (y0-ys)^2);
    numray{i} = numray1(i);
    if d0 <= r
        numray{i} = ray(i);
    elseif x0 < 0 && y0 > 0.2
        Rest = d0;
        if d0 <= 2*r
            Rest = 2*d0;
        end
        c0 = speed(node(i,:));
        Nray = 0;
        [ang] = NMLA_2D_2nd(x0,y0,c0,omega,Rest,rnode,relem,ru,ux,uy,pde,pct,Nray,data,opt,plt);
        if size(ang,2) > 2
            i
            ang
        end
        numray{i} = exp(1i*ang);
    end
end
toc;

figure(3)
ray_field(numray,node,6);



%% Ray_FEM
wpml = 6*(2*pi*cmin)/omega;        % width of PML
sigmaMax = 50/wpml; 
[u2] = Ray_FEM_PML_22_PointSource(node,elem,omega,wpml,sigmaMax,xs,ys,speed,numray,fquadorder,plt);
figure(4);
FJ_showresult(node,elem,real(u2));


[us] = Standard_FEM_PML_PointSource(node,elem,omega,wpml,sigmaMax,xs,ys,speed,fquadorder,plt);



%% Reference
fh = 1/2000;
[fnode,felem] = squaremesh([-0.5,0.5,-0.5,0.5],fh);
[uf] = Standard_FEM_PML_PointSource(fnode,felem,omega,wpml,sigmaMax,xs,ys,speed,fquadorder,plt);
figure(5);
FJ_showresult(fnode,felem,real(uf));

%% plot
a = 0.5;
[X,Y] = meshgrid(-a:h:a,-b:h:b);
[m,n] = size(X);

u2 = reshape(u2,m,n);
us = reshape(us,m,n);


fn = round(1/fh) + 1;
uf = reshape(uf,fn,fn);


figure(7);
hold off;
dy = 0.9;
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

figure(8);
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
plot(rxx,real(ruu),'k');
xlabel('x');
ylabel('Wavefield');
legend('Ray-FEM solution','Reference solution','LOCATION','Best');


figure(9);
hold off;
dy = 0.8;
yn = round(dy/h) + 1;
xx = X(yn,:);
uu = u2(yn,:);
ryn = round(dy/fh) + 1;
rxx = -0.5:fh:0.5;
ruu = uf(ryn,:);
ss = us(yn,:);
plot(xx,real(ss),'bs-');
hold on;
plot(rxx,real(ruu),'k');
xlabel('x');
ylabel('Wavefield');
legend('standard FEM solution','Reference solution','LOCATION','Best');


figure(10);
hold off;
dy = 0.9;
yn = round(dy/h) + 1;
xx = X(yn,:);
uu = u2(yn,:);
ryn = round(dy/fh) + 1;
rxx = -0.5:fh:0.5;
ruu = uf(ryn,:);
plot(xx,real(uu),'ro-');
hold on;
plot(rxx,real(ruu),'k');
xlabel('x');
ylabel('Wavefield');
legend('Ray-FEM solution','Reference solution','LOCATION','Best');


figure(11);
hold off;
dy = 0.9;
yn = round(dy/h) + 1;
xx = X(yn,:);
uu = u2(yn,:);
ryn = round(dy/fh) + 1;
rxx = -0.5:fh:0.5;
ruu = uf(ryn,:);
ss = us(yn,:);
plot(xx,real(ss),'bs-');
hold on;
plot(rxx,real(ruu),'k');
xlabel('x');
ylabel('Wavefield');
legend('standard FEM solution','Reference solution','LOCATION','Best');




end














%%
if (0);
    xs = -0.2;   ys = -0.3;             % point source location
speed = @(x) ( 1+ 0.5*sin(2*pi*x(:,1)));
cmin = 1/2;

Nray = 0;
Rest = 1;
pde = [];
pct = 1/2;
data = 'num';
opt = 0;

a= 1/2;
plt = 0;                           % plot solution or not
fquadorder = 3; 
omega = 60*pi;
high_omega = 60*pi;
minwavelength = (2*pi*cmin)/high_omega;

wpml = 6*minwavelength;                 % width of PML
sigmaMax = 25/wpml;                % Maximun absorbtion
r = 4*minwavelength;

NPW = 40;
rh = 1/(NPW*round(omega/(2*pi*cmin)));
fprintf('omega/(2*pi) = %d,   1/h = %d,   NPW = %d \n',omega/(2*pi), 1/rh, NPW);
[lnode,lelem] = squaremesh([-a,a,-a,a],rh);
load('test10_reference_solution_40.mat');
m = size(lnode,1);
m = round(sqrt(m));
n = m;

[ux,uy] = num_derivative_2(ru,rh,m,n,2);
a = 0.5;
b = 0.5;
mbd = 0;
% mbd = high_r;
% [mnode,melem] = squaremesh([-a-mbd,a+mbd,-b-mbd,b+mbd],h);
ch = 1/60;
[cnode,celem] = squaremesh([-a-mbd,a+mbd,-b-mbd,b+mbd],ch);
cN = size(cnode,1);
cnumray = cell(cN,1);
cray = ex_ray(cnode,xs,ys);

fprintf('NMLA time: \n');
tic;
for i = 1:cN
    x0 = cnode(i,1);  y0 = cnode(i,2);
    d0 = sqrt((x0-xs)^2 + (y0-ys)^2);
    if d0 < eps 
        cnumray{i} = 0;
    elseif d0 <= r
        cnumray{i} = exp(1i*cray(i,:));
    else
        Rest = d0;
        if d0 <= 2*r
            Rest = 2*d0;
        end
        c0 = speed(cnode(i,:));
        [ang] = NMLA_2D_2nd(x0,y0,c0,omega,Rest,lnode,lelem,ru,ux,uy,pde,pct,Nray,data,opt,plt);
        cnumray{i} = exp(1i*ang);
    end
end
toc;

ray_field(cnumray,cnode,2);
end

if (0)
%% Step 2: Use NMLA to find ray directions d_c with low frequency sqrt(\omega)
fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('\nStep2: NMLA, low frequency\n');

[ux,uy] = num_derivative_2(u_std,h,m,n,2);

mbd = high_r;
[mnode,melem] = squaremesh([-a-mbd,a+mbd,-b-mbd,b+mbd],h);
[cnode,celem] = squaremesh([-a-mbd,a+mbd,-b-mbd,b+mbd],ch);
cN = size(cnode,1);
cnumray = cell(cN,1);
cray = ex_ray(cnode,xs,ys);

fprintf('NMLA time: \n');
tic;
for i = 1:cN
    x0 = cnode(i,1);  y0 = cnode(i,2);
    d0 = sqrt((x0-xs)^2 + (y0-ys)^2);
    if d0 <= r
        cnumray{i} = exp(1i*cray(i,:));
    else
        Rest = d0;
        if d0 <= 2*r
            Rest = 2*d0;
        end
        c0 = speed(cnode(i,:));
        [ang] = NMLA_2D_2nd(x0,y0,c0,omega,Rest,lnode,lelem,u_std,ux,uy,pde,pct,Nray,data,opt,plt);
        cnumray{i} = exp(1i*ang);
    end
end
toc;

% clear lnode lelem;

% cdiffang = angle_error(cnumray,cray);
% norm(cdiffang,2)/norm(cray)
% norm(cdiffang,inf)

% cnumray = exp(1i*cnumray);

ray_field(cnumray,cnode,4);
% numray1 = interpolation(cnode,celem,mnode,cnumray);
% 
% ray = ex_ray(mnode,xs,ys);
% ray = exp(1i*ray);
% md = sqrt((mnode(:,1)-xs).^2 + (mnode(:,2)-ys).^2);
% ray = ray.*(1 - (md<eps));
% 
% numray1 = numray1.*(md>r) + ray.*(md<=r);
% 
% ray_field(numray1,node,10);


% diffray1 = numray1 - ray;
% figure(1);
% FJ_showresult(mnode,melem,real(diffray1));


end



%% ray direction field
if (0)

%% test for ray_fem_pml_2_pointsource
xs = 0;   ys = 0;  
speed = @(x) ones(size(x,1),1);    % medium speed

plt = 0;                           % plot solution or not
fquadorder = 3; 
omega = 20*pi;
NPW = 10;
h = 1/100;
% h = 1/(40*round(omega*NPW/(2*pi*40)));
wpml = 2*pi/omega;                 % width of PML
sigmaMax = 25/wpml;    

a = 1/2;
[node,elem] = squaremesh([-a,a,-a,a],h);
x = node(:,1);
y = node(:,2);

ray = ex_ray(node,xs,ys);
ray = exp(1i*ray);
d = sqrt((node(:,1)-xs).^2 + (node(:,2)-ys).^2);
ray = ray.*(1 - (d < eps) );
figure(1);
ray_field(ray,node,4);

% [u] = Ray_FEM_PML_1_PointSource(node,elem,omega,wpml,sigmaMax,xs,ys,speed,ray,fquadorder,plt);

N = size(node,1);
ray2 = cell(N,1);
for i = 1:N
    ray2{i} = ray(i);
end
% [u2,A2,v2,b2,f2] = Ray_FEM_PML_2_PointSource(node,elem,omega,wpml,sigmaMax,xs,ys,speed,ray2,fquadorder,plt);
figure(2);
ray_field(ray2,node,5);

% [u3] = Ray_FEM_PML_22_PointSource(node,elem,omega,wpml,sigmaMax,xs,ys,speed,ray2,fquadorder,plt);

end



if (0)
NPW = 40;
rh = 1/(NPW*round(omega/(2*pi*cmin)));
1/rh
[rnode,relem] = squaremesh([-a,a,-a,a],rh);
rN = size(rnode,1);
rn = round(sqrt(rN));
ru = reshape(u,rn,rn);

hold off;
dy = 0.9;
yn = round(dy/h) + 1;
xx = X(yn,:);
uu = uh(yn,:);
% suu = su(yn,:);
ryn = round(dy/rh) + 1;
rxx = -0.5:rh:0.5;
ruu = ru(ryn,:);
plot(xx,real(uu),'ro-');
hold on;
% plot(xx,real(suu),'b^-');
% hold on;
plot(rxx,real(ruu),'k');
xlabel('x');
ylabel('Wavefield');
legend('Ray-FEM solution','Standard FEM solution','Reference solution','LOCATION','Best');
% axis equal;

end


if (0)

xs = -1/8;   ys = -1/10;  
speed = @(x) ones(size(x,1),1);    % medium speed

plt = 0;                           % plot solution or not
fquadorder = 3; 
omega = 60*pi;
NPW = 20;
h = 1/(40*round(omega*NPW/(2*pi*40)));
wpml = 2*pi/omega;                 % width of PML
sigmaMax = 25/wpml;    

a = 1/2;
[node,elem] = squaremesh([-a,a,-a,a],h);
x = node(:,1);
y = node(:,2);

ray = ex_ray(node,xs,ys);
ray = exp(1i*ray);
d = sqrt((node(:,1)-xs).^2 + (node(:,2)-ys).^2);
ray = ray.*(1 - (d < eps) );

[u,~,v] = Ray_FEM_PML_1_PointSource(node,elem,omega,wpml,sigmaMax,xs,ys,speed,ray,fquadorder,plt);
% FJ_showresult(node,elem,real(u));
% showsolution(node,elem,real(u),2);



% fh = 1/1000;
% [fnode,felem] = squaremesh([-a,a,-a,a],fh);
% uray = ray_solution(node,elem,omega,speed,v,ray,fnode);
% figure(2);
% FJ_showresult(fnode,felem,real(uray));





r1 = 10.5*wpml;
r2 = 11.5*wpml;
theta1 = pi/8;
theta2 = 3*pi/8;

[theta,r] = meshgrid(theta1:1/10000:theta2, r1: 1/5000:r2);
[m,n] = size(r);
xx = r.*cos(theta) + xs;
yy = r.*sin(theta) + ys;
xynode = [xx(:),yy(:)];
uu = ray_solution(node,elem,omega,speed,v,ray,xynode);
uu = reshape(uu,m,n);
mesh(theta,r,real(uu))


end









% dist = sqrt((x-xs).^2 + (y-ys).^2);
% angle = atan2(y-ys, x-xs);
% angle = angle + 2*pi*(angle<0);
% index = (r1<=dist + eps).*(dist<=r2 + eps)...
%     .*(theta1 <= angle + eps).*(angle<=theta2 +eps);
% index = find(index);
% 
% theta = angle(index);
% r = dist(index);
% u_polar = real(u(index));
% 
% trisurf(delaunay(theta,r),theta,r,u_polar)
% 
% 






%% plot for vector field 
if (0)
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
quiver(cX,cY,cre,cim)
title('Ray direction field');


figure(3);
contour(X,Y,real(uh));
hold on
quiver(cX,cY,cre,cim)
title('Ray direction field');
end









%% test for point source
if (0)
global omega;
global a;
global h;


x0 = 0;   y0 = 0;
a = 1/2;
plt = 0;                   % show solution or not
fquadorder = 6;            % numerical quadrature order
omega = 20*pi;

wpml = 2*pi/omega;
sigmaMax = 25/wpml;

fh = 1/1000;
[fnode,felem] = squaremesh([-a,a,-a,a],fh);
speed = @(x) ones(size(x,1),1);



[u,A,b] = Standard_FEM_PML_PointSource(fnode,felem,omega,wpml,sigmaMax,x0,y0,speed,fquadorder,plt);

[X,Y] = meshgrid(-a:fh:a,-a:fh:a);
[m,n] = size(X);
uh = reshape(u,m,n);


figure(2)
contour(X,Y,real(uh));
title('level set of u_{std}');



ch = 1/200;
[node,elem] = squaremesh([-a,a,-a,a],ch);
p = node;
xx = p(:,1)-x0; yy = p(:,2)-y0;
ray = atan2(yy,xx);
ray = exp(1i*ray);
ray = ray.*(1 - (abs(xx)<eps).*(abs(yy)<eps));

[u2,A2,v2,b2] = Ray_FEM_PML_1_PointSource(node,elem,omega,wpml,sigmaMax,x0,y0,speed,ray,fquadorder,plt);


uh = u2;
[X,Y] = meshgrid(-a:ch:a,-a:ch:a);
[m,n] = size(X);
uh = reshape(uh,m,n);


figure(2)
contour(X,Y,real(uh));
title('level set of u_{ray}');



u1 = interpolation(fnode,felem,node,u);
norm(u1-u2,inf)
ch*norm(u2-u1,2)/(ch*norm(u1,2))


% x = node(:,1);  y = node(:,2);
% ue = 1i/4*besselh(0,1,omega*sqrt((x-x0).^2 + (y-y0).^2));
% norm(u - ue,inf)

end

if (0)
%% test7 NMLA for varying slowness
pde = Helmholtz_data7;
Nray = 1;
plt = 1;
global omega;
omega = 100*pi;
Rest = .35;
h = 1/400;
[node,elem] = squaremesh([-a,a,-a,a],h);

x0 = 0; y0 = 0;
c0 = pde.speed([x0,y0]);
[cnumray] = NMLA_2D_2nd(x0,y0,c0,omega,Rest,node,elem,0,0,0,pde,1/5,Nray,'ex',plt);
ray = pde.ray_ang([x0,y0]);
cnumray - ray
end



if (0)
%% test8 with exact ray
%% Load source data
pde = Helmholtz_data8;
xc = 1/8;   yc = 1/10;
Nray = 1;

%% Set up
plt = 0;                   % show solution or not
fquadorder = 6;            % numerical quadrature order
solver = 'DIR';  


a = 1/2;
omega = 20*pi;



NPW = 20;
wpml = 2*pi/omega;
sigmaMax = 25/wpml;

% fh = 1/(2^7*2^round(log2(omega/(2*pi))));
fh = 1/1600;
1/fh
[fnode,felem] = squaremesh([-a,a,-a,a],fh);
[u,~,~] = Standard_FEM_PML(fnode,felem,omega,wpml,sigmaMax,pde,fquadorder,solver,plt);



h = 1/(NPW*2^round(log2(omega/(2*pi))));
1/h
[snode,selem] = squaremesh([-a,a,-a,a],h);
sN = size(snode,1);
ray = pde.ray(snode);
[uh,~,v3,~] = Ray_FEM_PML_1(snode,selem,omega,wpml,sigmaMax,pde,ray,fquadorder,plt);

% ray = cell(sN,1);
% for i = 1:sN
%     xx = snode(i,1);  yy = snode(i,2);
%     if (abs(xx)<eps)*(abs(yy)<eps)
%         angl = linspace(0,2*pi,18) ;  
%         ang=angl(1:17) ;
%         ray{i} = exp(1i*ang);
%         i
%     else
%     ray{i} = pde.ray(snode);
%     end
% end
% 
% [uh,~,v3,~] = Ray_FEM_PML_2(snode,selem,omega,wpml,sigmaMax,pde,ray,fquadorder,plt);

ur = interpolation(fnode,felem,snode,u);

ue = pde.ex_u(snode);
diff1 = ue - uh;
diff2 = ur-uh;
h*norm(diff1,2)/(h*norm(ue,2))
h*norm(diff2,2)/(fh*norm(u,2))

[X,Y] = meshgrid(-a:h:a,-a:h:a);
[m,n] = size(X);
uh = reshape(uh,m,n);

figure(1);
FJ_showresult(snode,selem,real(diff2));
title('u_{std} - u_{ray}');

figure(2)
contour(X,Y,real(uh));
title('level set of u_{ray}');

figure(3);
contour(X,Y,real(reshape(ur,m,n)));
title('level set of u_{std}');

end

%% show plot for varing slowness case
if (0)
%     pde = Helmholtz_data7;
%     h = 1/1000;
%     a = 1/2;
%     [node,elem] = squaremesh([-a,a,-a,a],h);
%     ray_dir = pde.ray_ang(node);
%     angnorm = h*norm(ray_dir(:),2);
    omega = [ 10  20  30  40   ]*2*pi;
    ER_err = [ 6.72e-05  3.38e-05  2.26e-05  1.69e-05   ];
    NR_err = [ 1.05e-03  6.39e-04  3.16e-04  2.88e-04  ];
    ang_err = [4.79e-03  3.68e-03  3.25e-03  2.99e-03 ];
    
    figure(1);
    loglog(omega/(2*pi), ER_err,'bs-');
    hold on;
    loglog(omega/(2*pi), NR_err,'r*-');
    hold on;
    loglog(omega/(2*pi), ang_err,'ko-');
    %     hold on;
    %     loglog(omega/(2*pi), rec_int_err,'g^-');
    axis([0 40 -inf inf]);
    xlabel('Frequency \omega/2\pi');
    ylabel('Relative L^2 error');
    legend('ER-FEM Error','NR-FEM Error','Angle Error','LOCATION','Best');
    
    
    
    h = [1.00e-02  5.00e-03  2.50e-03  1.25e-03];
    NR_err = [1.13e-02  2.38e-03  6.04e-04  1.55e-04];
    ER_err = [5.22e-04  1.45e-04  3.86e-05  9.86e-06];
    
    figure(2);
    loglog(1./h, ER_err,'bs-');
    hold on;
    loglog(1./h, NR_err,'r*-');
    axis([100 1000 -inf inf]);
    xlabel('Mesh size 1/h');
    ylabel('Relative L^2 error');
    legend('ER-FEM Error','NR-FEM Error','LOCATION','Best');
    
end



%% Show convergence rate 
if (0)
    % load data
    load('rec1.mat');
    
    % angle error
    figure(1);
    loglog(rec_omega/(2*pi), rec_ang_err1,'bs-');
    hold on;
    loglog(rec_omega/(2*pi), rec_ang_err2,'r*-');
    axis([0 200 -inf inf])% axis tight;
    xlabel('Frequency \omega/2\pi');
    ylabel('Relative L^2 error');
    legend('Angle Error 1','Angle Error 2','LOCATION','Best');
    
    % NR-FEM error
    figure(2);
    loglog(rec_omega/(2*pi), rec_NR_err1,'bs-');
    hold on;
    loglog(rec_omega/(2*pi), rec_NR_err2,'r*-');
    hold on;
    loglog(rec_omega/(2*pi), rec_ER_err,'ko-');
    hold on;
    loglog(rec_omega/(2*pi), rec_int_err,'g^-');
    axis([0 200 -inf inf]);
    xlabel('Frequency \omega/2\pi');
    ylabel('Relative L^2 error');
    legend('NR-FEM Error 1','NR-FEM Error 2','ER-FEM Error','Interpoltation Error','LOCATION','Best');
    
    % optimality constant
    figure(3);
    plot(rec_omega/(2*pi), rec_NR_err2./rec_int_err,'*-');
    axis([0 200 -inf inf])% axis tight;
    xlabel('Frequency \omega/2\pi');
    ylabel('Optimality relation C');
end

%% show the convergence rate of second order ocrrection NMLA
if (0)
omega = [ 2.00e+01  2.20e+01  2.40e+01  2.60e+01  2.80e+01...
    3.00e+01  3.20e+01  3.40e+01  3.60e+01  3.80e+01...
    4.00e+01  4.20e+01  4.40e+01  4.60e+01  4.80e+01...
    5.00e+01  5.20e+01  5.40e+01  5.60e+01  5.80e+01...
    6.00e+01  6.20e+01  6.40e+01  6.60e+01  6.80e+01...
    7.00e+01  7.20e+01  7.40e+01  7.60e+01  7.80e+01...
    8.00e+01  ]*2*pi;

ang_err = [ 1.28e-05  1.30e-05  6.41e-06  2.11e-05  1.14e-05...
    8.10e-06  1.63e-05  9.25e-06  1.18e-05  5.84e-06...
    5.16e-06  8.86e-06  4.54e-06  4.42e-06  5.75e-06...
    3.67e-06  4.01e-06  3.57e-06  4.76e-06  3.71e-06...
    4.39e-06  3.43e-06  3.37e-06  2.79e-06  2.38e-06...
    2.78e-06  3.83e-06  3.12e-06  3.18e-06  2.31e-06...
    3.79e-06];

NR_err = [1.96e-05  1.71e-05  1.63e-05  1.53e-05  1.43e-05...
    1.33e-05  1.25e-05  1.19e-05  1.13e-05  1.07e-05...
    1.00e-05  9.53e-06  9.09e-06  8.76e-06  8.36e-06...
    8.11e-06  7.75e-06  7.54e-06  7.24e-06  7.03e-06...
    6.67e-06  6.50e-06  6.14e-06  5.84e-06  5.73e-06...
    5.65e-06  5.53e-06  5.41e-06  5.30e-06  5.17e-06...
    5.06e-06  ];

ER_err = [2.00e-05  1.82e-05  1.67e-05  1.54e-05  1.43e-05...
    1.33e-05  1.25e-05  1.18e-05  1.11e-05  1.05e-05...
    1.00e-05  9.54e-06  9.11e-06  8.71e-06  8.35e-06...
    8.02e-06  7.71e-06  7.43e-06  7.16e-06  6.92e-06...
    6.69e-06  6.47e-06  6.27e-06  6.08e-06  5.90e-06...
    5.73e-06  5.57e-06  5.42e-06  5.28e-06  5.14e-06...
    5.02e-06 ];


figure(1);
FJ_showrate(omega,ang_err)
figure(2);
FJ_showrate(omega,NR_err)
figure(3);
semilogy(omega/(2*pi),ang_err,'r*-');
hold on;
semilogy(omega/(2*pi),NR_err,'s-');
end


%% test for the real radius
if (0)
    Rest = 2;
    n = 10;
    omega = zeros(n,1);
    r = zeros(n,1);
    for i = 1:n
        omega(i) = (i+4)*2*pi;
        p = [1,0,1,-1.5-0.59*(omega(i)*Rest)^0.75];
        rt = roots(p);
        pidx = find(rt>0);
        r(i) = (rt(pidx(1)))^3/omega(i);
    end
    
    showrate(omega,r)
end




%% test for second order correction
if (0)
    pde = Helmholtz_data1;
    x0 = 15;
    y0 = 1.1;
    omega = 320*pi;
    Rest = sqrt((x0-2)^2 + (y0-2)^2);
    % [numray] = NMLA_2D(x0,y0,1,omega,Rest,0,0,0,0,0,pde,1/5,1,'ex',1);
    [numray] = NMLA_2D_2nd(x0,y0,1,omega,Rest,0,0,0,0,0,pde,1/5,1,'ex',1);
    
    ray = pde.ray([x0,y0]);
    ray = [real(ray), imag(ray)];
    ray_dir = atan2(ray(:,2),ray(:,1));
    ray_dir = ray_dir + 2*pi*(ray_dir<0);
    
    numray - ray_dir
    
    
    xx = 1:10;
    yy = 2*xx.^2 + 4*xx + 3;
    p = polyfit(xx,yy,2);
end


%% show the complexity of interpolation
if (0)
    n = 6;
    hi = 1:n;
    h = 1./(2.^(hi+8));
    h1 = 1./(2.^(hi+2));
    N = zeros(n,1);
    NT = N;
    
    for i = 1:n
        [node,elem] = squaremesh([0,1,0,1],h(3));
        [cnode,celem] = squaremesh([0,1,0,1],h1(i));
        cN = size(cnode,1);
        u = rand(cN,1);
        tic;
        uint = interpolation(cnode,celem,node,u);
        time = toc
        size(node,1)
    end
end

%% show the rate of NMLA radius wrt frequency
if (0)
    Rest = 4;
    n = 1000;
    
    omega = 1:n;
    omega = omega'*2*pi;
    pr = zeros(n,1);
    for ii = 1: n
        p = [1,1,-2.5-(omega(ii)*Rest)];
        rt = roots(p);
        pidx = find(rt>0);
        pr(ii) = rt(pidx(1))/omega(ii);
    end
    showrate(omega,pr)
end
