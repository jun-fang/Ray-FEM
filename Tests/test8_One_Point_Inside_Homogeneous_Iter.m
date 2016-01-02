%% One point source inside domain (homogeneous medium):  iterative idea

xs = -0.3;   ys = -0.3;             % point source location
speed = @(x) ones(size(x,1),1);    % medium speed

plt = 0;                           % plot solution or not
fquadorder = 3;                    % numerical quadrature order
Nray = 1;
Rest = 1;
pde = [];
pct = 1/5;
data = 'num';
opt = 1;

high_omega = 200*pi;
low_omega = sqrt(high_omega);
NPW = 10;

h = 1/(20*round(high_omega*NPW/(2*pi*20)));
ch = 1/(20*max(round(low_omega*NPW/(2*pi*20)),1));

wavelength = 2*pi/high_omega;  
wpml = wavelength;               % width of PML
sigmaMax = 25/wpml;                % Maximun absorbtion
r = 8*wpml;

lg_a = 1;
md_a = 0.65;

sm_a = 1/2;

fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('omega/(2*pi) = %d,   1/h = %d   1/ch = %d,  NPW = %d \n',high_omega/(2*pi), 1/h, 1/ch, NPW);


%% Step 1: Solve the Hemholtz equation with the same source but with a relative low frequency sqrt(\omega) by Standard FEM, mesh size \omega*h = constant
fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('Step1: S-FEM, low frequency\n');

a = lg_a;
[lnode,lelem] = squaremesh([-a,a,-a,a],h);
omega = low_omega;
[u_std] = Standard_FEM_PML_PointSource(lnode,lelem,omega,wpml,sigmaMax,xs,ys,speed,fquadorder,plt);


%% Step 2: Use NMLA to find ray directions d_c with low frequency sqrt(\omega)
fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('\nStep2: NMLA, low frequency\n');

[ux,uy] = num_derivative(u_std,h,2);

a = md_a;
[mnode,melem] = squaremesh([-a,a,-a,a],h);
[cnode,celem] = squaremesh([-a,a,-a,a],ch);
cN = size(cnode,1);
cnumray = zeros(cN,Nray);
cray = ex_ray(cnode,xs,ys);

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

cdiffang = angle_error(cnumray,cray);
norm(cdiffang,2)/norm(cray)
norm(cdiffang,inf)

cnumray = exp(1i*cnumray);
numray1 = interpolation(cnode,celem,mnode,cnumray);

ray = ex_ray(mnode,xs,ys);
ray = exp(1i*ray);
md = sqrt((mnode(:,1)-xs).^2 + (mnode(:,2)-ys).^2);
ray = ray.*(1 - (md<eps));

numray1 = numray1.*(md>r) + ray.*(md<=r);
diffray1 = numray1 - ray;
% figure(1);
% FJ_showresult(mnode,melem,real(diffray1));

%% Step 3: Solve the original Helmholtz equation by Ray-based FEM with ray directions d_c
fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('\nStep3: Ray-FEM, high frequency \n');

omega = high_omega;
[u1] = Ray_FEM_PML_1_PointSource(mnode,melem,omega,wpml,sigmaMax,xs,ys,speed,numray1,fquadorder,plt);


%% Step 4: NMLA to find original ray directions d_o with wavenumber k
fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('\nStep4: NMLA, high frequency \n');

a = sm_a;
[node,elem] = squaremesh([-a,a,-a,a],h);
[cnode,celem] = squaremesh([-a,a,-a,a],ch);
cN = size(cnode,1);
cnumray = zeros(cN,Nray);
cray = ex_ray(cnode,xs,ys);

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

cdiffang = angle_error(cnumray,cray);
norm(cdiffang,2)/norm(cray)
norm(cdiffang,inf)

cnumray = exp(1i*cnumray);
numray2 = interpolation(cnode,celem,node,cnumray);

ray = ex_ray(node,xs,ys);
ray = exp(1i*ray);
d = sqrt((node(:,1)-xs).^2 + (node(:,2)-ys).^2);
ray = ray.*(1 - (d < eps) );

numray2 = numray2.*(d>r) + ray.*(d<=r);
diffray2 = numray2 - ray;
% figure(2);
% FJ_showresult(node,elem,real(diffray2));

numray_dir = [real(numray2), imag(numray2)];
numray_dir = atan2(numray_dir(:,2), numray_dir(:,1));
numray_dir = numray_dir + 2*pi*(numray_dir<0);

diffang = angle_error(numray_dir,ex_ray(node,xs,ys));
diffang = diffang.*(d>r);
% figure(3);
% FJ_showresult(node,elem,real(diffang));
% title('NMLA angle error');


%% Step 5: Solve the original Helmholtz equation by Ray-based FEM with ray directions d_o
fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('\nStep5: Ray-FEM, high frequency \n');

omega = high_omega;
[u2,~,v2] = Ray_FEM_PML_1_PointSource(node,elem,omega,wpml,sigmaMax,xs,ys,speed,numray2,fquadorder,plt);


%% 
% [X,Y] = meshgrid(-a:h:a,-a:h:a);
% [m,n] = size(X);
% uh = reshape(u2,m,n);
% 
% figure(4)
% contour(X,Y,real(uh));
% title('Level set of Ray-FEM solution u_{ray}');


%% map to polar
figure(5);
% r1 = 10.3*wpml;
% r2 = 12.5*wpml;
r1 = 0.834;
r2 = 0.834 + 2.1*wavelength;
theta1 = pi/4 - pi/16;
theta2 = pi/4 + pi/16;
subplot(2,1,1);
mapto_polar(node,elem,omega,speed,v2,numray2,xs,ys,r1,r2,theta1,theta2);

subplot(2,1,2);
mapto_polar(node,elem,omega,speed,v2,numray2,xs,ys,r1,r2,theta1,theta2);
axis equal



%% 
% figure(6);
% % r1 = 0.81;
% % r2 = 0.81 + 2.1*wavelength;
% % omega = 80*pi;
% [X,Y] = meshgrid(-a:1/1000:a,-a:1/1000:a);
% [m,n] = size(X);
% xynode = [X(:),Y(:)];
% uu = ray_solution(node,elem,omega,speed,v2,numray2,xynode);
% uu = reshape(uu,m,n);
% mesh(X,Y,real(uu));
% az = 0;
% el = 90;
% view(az, el);
% axis equal; axis tight;
% hold on;
% 
% theta = theta1:1/10000:theta2;
% rr = r1: 1/5000:r2;
% x1 = r1*cos(theta) + xs;   y1 = r1*sin(theta) + ys;
% x2 = r2*cos(theta) + xs;   y2 = r2*sin(theta) + ys;
% x3 = rr*cos(theta1) + xs;  y3 = rr*sin(theta1) + ys;
% x4 = rr*cos(theta2) + xs;  y4 = rr*sin(theta2) + ys;
% p = plot(x1,y1,'r-');
% p.LineWidth = 3;  hold on;
% p = plot(x2,y2,'r-');
% p.LineWidth = 3;  hold on;
% p = plot(x3,y3,'r-');
% p.LineWidth = 3;  hold on;
% p = plot(x4,y4,'r-');
% p.LineWidth = 3;  hold on;
% xlabel('x');
% ylabel('y');
