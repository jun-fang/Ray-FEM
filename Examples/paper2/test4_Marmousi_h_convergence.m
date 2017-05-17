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

% domain
sdx = 1.5; sdy = 0.5;
small_domain = [-sdx, sdx, -sdy, sdy];

% real frequency in Hertz
real_omega = high_omega*1500/4000;
hz = real_omega/(2*pi);

% print
fprintf(['-'*ones(1,80) '\n']);
fprintf('Marmousi wave speed case at %.2f Hz: \n', hz);
fprintf(['-'*ones(1,80) '\n']);
fprintf('Computational domain = \n  [%.2f, %.2f, %.2f, %.2f] \n', small_domain);


%% load reference data
load('results_4_Marmousi_100pi_NPW_16.mat');
ray_ref = ray;  uh_ref = uh2;  h_ref = 1/800;
[rnode, relem] = squaremesh(small_domain, h_ref);
rx = rnode(:,1); ry = rnode(:,2);
rr = sqrt((rx-xs).^2 + (ry-ys).^2);
ruu = uh_ref;  ruu(rr <= 2*wl) = 0;
rnorm = norm(ruu)*h_ref;

% mesh size
h = 1/25;

test_num = 3;
err = zeros(test_num,1);
hs = err;

for j = 1:test_num
    
    h = h/2;
    hs(j) = h;
    tic;
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
    toc;
    
    fprintf(['-'*ones(1,80) '\n']);
    fprintf('  Wavelength = %.2d    NPW = %d    1/h = %d  \n', wl, round(wl/h), round(1/h) );
    
    
    
    tstart = tic;
    [node,elem] = squaremesh(small_domain,h);
    N = size(node,1);  Ndof = 0;
    ray = cell(N,1);  ur = zeros(N,1);
    
    nh = h/h_ref;
    
    m_ref = round( (small_domain(2) - small_domain(1)) /h_ref ) + 1;
    n_ref = round( (small_domain(4) - small_domain(3)) /h_ref ) + 1;
    
    m = round( (small_domain(2) - small_domain(1)) /h ) + 1;
    n = round( (small_domain(4) - small_domain(3)) /h ) + 1;
    
    tic;
    for i = 1:N
        mi = ceil(i/n);
        ni = i - (mi-1)*n;
        mi_ref = (mi-1)*nh + 1;
        ni_ref = (ni-1)*nh + 1;
        i_ref = (mi_ref-1)*n_ref + ni_ref;
        ray{i} = ray_ref{i_ref};
        Ndof = Ndof + length(ray{i});
        ur(i) = uh_ref(i_ref);
    end
    toc;
    
    % figure(1);
    % subplot(2,1,1);
    % showsolution(node,elem,real(uh),2);
    %
    % [rnode,relem] = squaremesh(small_domain,h_ref);
    %
    % subplot(2,1,2);
    % showsolution(rnode,relem,real(uh_ref),2);
    
    % figure(74); ray_field(ray,node,20,1/10); axis equal; axis tight;
    
    
    %% Step 5: Solve the original Helmholtz equation by Ray-based FEM with ray directions d_o
    fprintf([ '-'*ones(1,80) '\n']);
    fprintf('Step5: Ray-FEM, high frequency \n');
    
    omega = high_omega;
    wpml = 0.1;                % width of PML
    sigmaMax = 25/wpml;                 % Maximun absorbtion
    option ='homogeneous';
    
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
    
    du = ur - uh2;
    du(rr <= 2*wl) = 0;
    
    err(j) = norm(du)*h/rnorm
    
end


totaltime = toc(tstart);
fprintf('\n\nTotal running time: % d minutes \n', totaltime/60);



showrate(hs,err);



% nameFile = strcat('results_4_Marmousi_',num2str(round(omega/pi)), 'pi_NPW_',num2str(round(wl/h)),'.mat');
% save(nameFile, 'uh2', 'h', 'high_omega','ray');
%
% % figure(70);
% % m = round( (small_domain(2) - small_domain(1)) /h ) + 1;
% % n = round( (small_domain(4) - small_domain(3)) /h ) + 1;
% % uh22 = reshape(real(uh2), n, m);
% % imagesc(uh22(end:-1:1, :)); axis equal; axis tight; colorbar;
%
% figure(75);
% figure('position', [100, 100, 1600, 525]);
% showsolution(node,elem,real(uh2),2); colorbar; axis equal; axis tight;
% % set(gca,'position',[0.05 0.05 .875 .95],'units','normalized')
% % set(gca,'position',[0.025 0.025 .915 0.875],'units','normalized')
% set(gca,'position',[0.01 0.075 .95 0.875],'units','normalized')


