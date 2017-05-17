%% Marmousi wave speed case

clear;
addpath(genpath('../../../ifem/'));
addpath('../../Functions/')
addpath('../../Methods/');
addpath('../../NMLA/');
addpath('../../Plots/');
addpath('../../Helmholtz_data/');
addpath('C:\Users\Jun Fang\Documents\MATLAB\Marmousi\');


%% Set up
plt = 0;                   % show solution or not
fquadorder = 3;            % numerical quadrature order
Nray = 4;                  % one ray direction
sec_opt = 0;               % NMLA second order correction or not
pct = 0.25;

xs = 0; ys = 0.3;          % source location
epsilon = 1/(4*pi);        % cut-off parameter

% frequency
high_omega = 100*pi;
low_omega = 2*sqrt(high_omega);
wl = 2*pi/high_omega;


% width of PML
high_wpml = 0.1;
low_wpml = 0.25;

% mesh size
h = 1/800;  

% load Marmousi data
load('Marmousi_smoother.mat');  
hr = 1/4000; mr = 8001; nr = 16001; xr = 2; yr = 1;

% compress Marmousi data
nh = round(h/hr);
ix = 1:nh:mr;  iy = 1:nh:nr;
Marmousi_compressed = Marmousi_smoother(ix,iy);
Marmousi_speed = Marmousi_compressed(:);
clear ix iy Marmousi_smoother Marmousi_compressed;

% construct Marmousi speed
speed = @(p) Marmousi_speed( Marmousi_index(p, xr, yr, h) )/1500;    % wave speed

% domain
sdx = 1.5; sdy = 0.5;
mdx = 1.65; mdy = 0.65;
ldx = 2; ldy = 1;

large_domain = [-ldx, ldx, -ldy, ldy];
middle_domain = [-mdx, mdx, -mdy, mdy];
small_domain = [-sdx, sdx, -sdy, sdy];


% real frequency in Hertz
real_omega = high_omega*1500/4000;
hz = real_omega/(2*pi);

% print 
fprintf(['-'*ones(1,80) '\n']);
fprintf('Marmousi wave speed case at %.2f Hz: \n', hz);
fprintf(['-'*ones(1,80) '\n']);
fprintf('Computational domain = \n  [%.2f, %.2f, %.2f, %.2f] \n', large_domain);
fprintf(['-'*ones(1,80) '\n']);
fprintf('  Wavelength = %.2d    NPW = %d    1/h = %d  \n', wl, round(wl/h), round(1/h) );


tstart = tic;
%% Step 1: Solve the Hemholtz equation with the same source but with a relative low frequency sqrt(\omega) by Standard FEM, mesh size \omega*h = constant
fprintf(['-'*ones(1,80) '\n']);
fprintf('Step1: S-FEM, low frequency \n');

tic;
omega = low_omega;              % low frequency
wpml = low_wpml;                % width of PML
sigmaMax = 25/wpml;                 % Maximun absorbtion
[lnode,lelem] = squaremesh(large_domain,h);

% smooth part
A = assemble_Helmholtz_matrix_SFEM(lnode,lelem,omega,wpml,sigmaMax,speed,fquadorder);
b = assemble_RHS_SFEM_with_ST(lnode,lelem,xs,ys,omega,wpml,sigmaMax,epsilon,fquadorder);
[~,~,isBdNode] = findboundary(lelem);
freeNode = find(~isBdNode);
lN = size(lnode,1);        u_std = zeros(lN,1);
u_std(freeNode) = A(freeNode,freeNode)\b(freeNode);

% singular part
x = lnode(:,1); y = lnode(:,2);
rr = sqrt((x-xs).^2 + (y-ys).^2);
ub = 1i/4*besselh(0,1,omega*rr);
cf = cutoff(epsilon,2*epsilon,lnode,xs,ys);

% low frequency solution: smooth + singularity
u_low = u_std + ub.*cf;
toc;

% figure(71); showsolution(lnode,lelem,real(u_low),2); colorbar;axis equal; axis tight;

clear A b cf freeNode isBdNode rr u_std ub x y;

%% Step 2: Use NMLA to find ray directions d_c with low frequency sqrt(\omega)
fprintf(['-'*ones(1,80) '\n']);
fprintf('Step2: NMLA, low frequency \n');

% compute numerical derivatives
m = round( (large_domain(2) - large_domain(1)) /h ) + 1;
n = round( (large_domain(4) - large_domain(3)) /h ) + 1;
[ux,uy] = num_derivative(u_low,h,2,m,n);

% all middle domain: homogeneous + inhomogeneous
[mnode,melem] = squaremesh(middle_domain,h);
mN = size(mnode,1);  mNdof = 0;  mcompressed = 0;
mray = cell(mN,1);

tic;
for mi = 1:mN
    x0 = mnode(mi,1);  y0 = mnode(mi,2);
    r0 = sqrt((x0-xs)^2 + (y0-ys)^2);
    c0 = speed(mnode(mi,:));
    Rest = min(1.45, r0);
    if r0 < 0.1  % near source 
        mray{mi} = ex_ray([x0,y0],xs,ys,1);
        mNdof = mNdof + 1;
    else  % far field
        angles = NMLA(x0,y0,c0,omega,Rest,lnode,lelem,u_low,ux,uy,[],1/4,Nray,'num',sec_opt,plt);
        mray{mi} = exp(1i*angles);
        mNdof = mNdof + length(angles);  
    end
end
toc;

% figure(72); ray_field(mray,mnode,20,1/10);
    
clear lnode lelem u_low ux uy;

%% Step 3: Solve the original Helmholtz equation by Ray-based FEM with ray directions d_c
fprintf(['-'*ones(1,80) '\n']);
fprintf('Step3: Ray-FEM, high frequency \n');

omega = high_omega;
wpml = high_wpml;                % width of PML
sigmaMax = 25/wpml;                 % Maximun absorbtion
ray = mray;

% smooth part
option ='homogeneous'; 
tic;
A = assemble_Helmholtz_matrix_RayFEM(mnode,melem,omega,wpml,sigmaMax,speed,ray,fquadorder);
toc;
tic;
b = assemble_RHS_RayFEM_with_ST(mnode,melem,xs,ys,omega,epsilon,wpml,sigmaMax,ray,speed,fquadorder,option);
toc;
tic;
uh = RayFEM_direct_solver(mnode,melem,A,b,omega,ray,speed);
toc;

% singularity part
x = mnode(:,1); y = mnode(:,2);
rr = sqrt((x-xs).^2 + (y-ys).^2);
ub = 1i/4*besselh(0,1,omega*rr);
cf = cutoff(epsilon,2*epsilon,mnode,xs,ys);

% smooth + singularity
uh1 = uh + ub.*cf;

% figure(73); showsolution(mnode,melem,real(uh1),2); colorbar; axis equal; axis tight;

clear A b cf mray rr ub uh x y; 

%% Step 4: NMLA to find original ray directions d_o with wavenumber k
fprintf(['-'*ones(1,80) '\n']);
fprintf('Step4: NMLA, high frequency \n');
tic;

% compute numerical derivatives
m = round( (middle_domain(2) - middle_domain(1)) /h ) + 1;
n = round( (middle_domain(4) - middle_domain(3)) /h ) + 1;
[ux,uy] = num_derivative(uh1,h,2,m,n);

% all small domain: homogeneous + inhomogeneous
[node,elem] = squaremesh(small_domain,h);
N = size(node,1);  Ndof = 0;  compressed = 0;
ray = cell(N,1);

tic;
for i = 1:N
    x0 = node(i,1);  y0 = node(i,2);
    r0 = sqrt((x0-xs)^2 + (y0-ys)^2);
    c0 = speed(node(i,:));
    Rest = min(1.5, r0);
    if r0 < 0.1  % near source 
        ray{i} = ex_ray([x0,y0],xs,ys,1);
        Ndof = Ndof + 1;
    else  % far field
        angles = NMLA(x0,y0,c0,omega,Rest,mnode,melem,uh1,ux,uy,[],1/4,Nray,'num',sec_opt,plt);
        ray{i} = exp(1i*angles);
        Ndof = Ndof + length(angles);  
    end
end
toc;

% figure(74); ray_field(ray,node,20,1/10); axis equal; axis tight;

clear mnode melem uh1 ux uy;

%% Step 5: Solve the original Helmholtz equation by Ray-based FEM with ray directions d_o
fprintf([ '-'*ones(1,80) '\n']);
fprintf('Step5: Ray-FEM, high frequency \n');

omega = high_omega;
wpml = 0.1;                % width of PML
sigmaMax = 25/wpml;                 % Maximun absorbtion

% Assembling
tic;
A = assemble_Helmholtz_matrix_RayFEM(node,elem,omega,wpml,sigmaMax,speed,ray,fquadorder);
toc; 
tic;
b = assemble_RHS_RayFEM_with_ST(node,elem,xs,ys,omega,epsilon,wpml,sigmaMax,ray,speed,fquadorder,option);
toc;
tic;
uh = RayFEM_direct_solver(node,elem,A,b,omega,ray,speed);
toc;

% singularity part
x = node(:,1); y = node(:,2);
rr = sqrt((x-xs).^2 + (y-ys).^2);
ub = 1i/4*besselh(0,1,omega*rr);
cf = cutoff(epsilon,2*epsilon,node,xs,ys);

% smooth + singularity
uh2 = uh + ub.*cf;

% figure(75); showsolution(node,elem,real(uh2),2); colorbar; axis equal; axis tight;


totaltime = toc(tstart);
fprintf('\n\nTotal running time: % d minutes \n', totaltime/60);

nameFile = strcat('results_4_Marmousi_',num2str(round(omega/pi)), 'pi_NPW_',num2str(round(wl/h)),'.mat');
save(nameFile, 'uh2', 'h', 'high_omega','ray','mNdof', 'Ndof');

% figure(70);
% m = round( (small_domain(2) - small_domain(1)) /h ) + 1;
% n = round( (small_domain(4) - small_domain(3)) /h ) + 1;
% uh22 = reshape(real(uh2), n, m);
% imagesc(uh22(end:-1:1, :)); axis equal; axis tight; colorbar;

% figure(75); 
% figure('position', [100, 100, 1600, 525]);
% showsolution(node,elem,real(uh2),2); colorbar; axis equal; axis tight;
% % set(gca,'position',[0.05 0.05 .875 .95],'units','normalized')
% % set(gca,'position',[0.025 0.025 .915 0.875],'units','normalized')
% set(gca,'position',[0.01 0.075 .95 0.875],'units','normalized')


