%% One point source inside domain: homogeneous medium

% point source information
xs = -0.4;   ys = -0.4;             % point source location
speed = @(x) ones(size(x,1),1);     % medium speed

% parameters for numerical methods
plt = 0;                            % plot solution or not
fquadorder = 6;                     % numerical quadrature order
Nray = 1;
Rest = 1;
pde = [];
pct = 1/5;
data = 'num';
opt = 1;                            % NMLA second order correction

% physical parameters
high_omega = 80*pi;                % high frequency
low_omega = sqrt(high_omega);       % low frequency
NPW = 10;                            % number of grid points per wavelength

starttime = cputime;
h = 1/(10*round(high_omega*NPW/(2*pi*10)));             % fine mesh
ch = 1/(10*max(round(low_omega*NPW/(2*pi*10)),1));      % coarse mesh

wl = 2*pi/high_omega;               % wavelength
wpml = h*ceil(wl/h);                % width of PML
sigmaMax = 25/wpml;                 % Maximun absorbtion
r = 8*wpml;

lg_a = 1.1;                           % large domain
md_a = 0.7;                        % middle domain
sm_a = 1/2;                         % small doamin

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
cray = ex_ray_angle(cnode,xs,ys);
% cray = exp(1i*cray);
cr1 = zeros(cN,1);
cd1 = cr1;

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
        [cnumray(i,:),cr1(i)] = NMLA_2D_2nd(x0,y0,c0,omega,Rest,lnode,lelem,u_std,ux,uy,pde,pct,Nray,data,opt,plt);
        cd1(i) = cr1(i) - d0;
        %         if r1 > d0
        %             cr1(i) = 1;
        %         end
    end
end
toc;
  
clear lnode lelem;

cdiffang = angle_error(cnumray,cray);
norm(cdiffang,2)/norm(cray)
norm(cdiffang,inf)

cnumray = exp(1i*cnumray);
numray1 = interpolation(cnode,celem,mnode,cnumray);

ray = ex_ray_angle(mnode,xs,ys);
ray = exp(1i*ray);
md = sqrt((mnode(:,1)-xs).^2 + (mnode(:,2)-ys).^2);
ray = ray.*(1 - (md<eps));

numray1 = numray1.*(md>r) + ray.*(md<=r);
diffray1 = numray1 - ray;
% figure(1);
% FJ_showresult(mnode,melem,real(diffray1));
% figure(1);
% FJ_showresult(cnode,celem,cr1);

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
cray = ex_ray_angle(cnode,xs,ys);
cr2 = zeros(cN,1);
cd2 = cr2;

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
        [cnumray(i,:),cr2(i)] = NMLA_2D_2nd(x0,y0,c0,omega,Rest,mnode,melem,u1,ux,uy,pde,pct,Nray,data,opt,plt);
        cd2(i) = cr2(i) - d0;
        %         if cr > d0
        %             cr2(i) = 1;
        %         end
    end
end
toc;

% clear mnode melem;

cdiffang = angle_error(cnumray,cray);
norm(cdiffang,2)/norm(cray)
norm(cdiffang,inf)

cnumray = exp(1i*cnumray);
numray2 = interpolation(cnode,celem,node,cnumray);

ray = ex_ray_angle(node,xs,ys);
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

diffang = angle_error(numray_dir,ex_ray_angle(node,xs,ys));
diffang = diffang.*(d>r);
% figure(3);
% FJ_showresult(node,elem,real(diffang));
% title('NMLA angle error');
% figure(3);
% FJ_showresult(cnode,celem,cr2);

%% Step 5: Solve the original Helmholtz equation by Ray-based FEM with ray directions d_o
fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('\nStep5: Ray-FEM, high frequency \n');

omega = high_omega;
[u,~,v2] = Ray_FEM_PML_1_PointSource(node,elem,omega,wpml,sigmaMax,xs,ys,speed,numray2,fquadorder,plt);
% [u] = Ray_FEM_PML_1_PointSource(node,elem,omega,wpml,sigmaMax,xs,ys,speed,ray,fquadorder,plt);


%% Exact Ray-FEM
% [ue,~,ve] = Ray_FEM_PML_1_PointSource(node,elem,omega,wpml,sigmaMax,xs,ys,speed,ray,fquadorder,plt);

%% S-FEM
% [us] = Standard_FEM_PML_PointSource(node,elem,omega,wpml,sigmaMax,xs,ys,speed,fquadorder,plt);



endtime = cputime;

fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('\nTotal running time: %d seconds\n', endtime-starttime);



%%
figure(10);
ray_field(numray2,node,10,1/5000);
axis([-0.6,0.6,-0.6,0.6]);
set(gca, 'FontSize', 18); 
set(gcf, 'Position', [250 250 550 500])
% axis equal; %axis tight;


%%
figure(20);
hold off;
% sp1 = subplot(2,1,1);
positionVector1 = [0.175, 0.275, 0.7, 0.7];
subplot('Position',positionVector1)

showsolution(node,elem,real(u),2)
axis equal; axis off; 
colorbar;
set(gca, 'FontSize', 18); 

hold on

nr = 30;
r1 = wl*(nr+0.4);
r2 = wl*(nr+2.4);
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

positionVector2 = [0.15, 0.15, 0.7, 0.1];
subplot('Position',positionVector2)

% sp2 = subplot(2,1,2);
polar_plot_solution(omega,node,xs,ys,r1,r2,theta1,theta2,u,'spline',1);
% axis equal;
axis tight;
% colorbar;
set(gca, 'FontSize', 18); 
set(gcf, 'Position', [1000 250 600 500])

