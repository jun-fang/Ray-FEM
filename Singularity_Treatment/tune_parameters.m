%% Tune parameters: 
% Domains: frequency --> PML --> ld, md

% NMLA: Rest --> NMLA radius --> ld, md
%       h_c, M --> convergence rate


%% Basic set up

% add paths
clear;
addpath(genpath('../../ifem/'));
addpath('../Methods/');
addpath('../NMLA/');
addpath('../Cutoff_Functions/')
addpath('../Plots_Prints/');

% load source data
xs = 0; ys = 0;                      % point source location
epsilon = 50/(80*pi);                % cut-off parameter   
speed = @(x) ones(size(x,1),1);      % medium speed

% high frequency 
high_omega = [120 160 200 240 320 480 640 800 960]*pi;

% others
plt = 0;                   % show solution or not
fquadorder = 3;            % numerical quadrature order
Nray = 1;                  % one ray direction at each grid node
sec_opt = 0;               % NMLA second order correction or not

% error
low_max_rayerr = 0*high_omega;            % L_inf ray error of low-freq waves
low_l2_rayerr = 0*high_omega;             % L_2 ray error of low-freq waves
high_max_rayerr = 0*high_omega;           % L_inf ray error of high-freq waves
high_l2_rayerr = 0*high_omega;            % L_2 ray error of high-freq waves

sfem_max_err = 0*high_omega;              % SFEM L_inf error of low-freq solution 
sfem_rel_max_err = 0*high_omega;          % SFEM relative L_inf error 
sfem_l2_err = 0*high_omega;               % SFEM L_2 error 
sfem_rel_l2_err = 0*high_omega;           % SFEM relative L_2 error 

rayfem_max_err = 0*high_omega;            % L_inf error of the numerical solution
rayfem_rel_max_err = 0*high_omega;        % relative L_inf error of the numerical solution
rayfem_l2_err = 0*high_omega;             % L_2 error of the numerical solution
rayfem_rel_l2_err = 0*high_omega;         % relative L_2 error of the numerical solution


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The choice of low frequency, PML
% L2 error should be O(omega^-3/4), relative L2 error should be O(\omega^-1/2)
% choose low_omega and low_pml such that the domain is as small as
% posssible
% low_omega = 2*sqrt(high_omega); low_wpml = 0.1 can achieve a small domain
% when NPW = 4, 6 or 8 fixed, the length of PML doesn't fact too much: 
% varying low_wpml = 0.1, 0.2, 0.3, the results are almost the same
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

low_omega = 2*sqrt(high_omega); 
NPW = 8;
test_num = 4;

fh = 1./(NPW*round(high_omega/(2*pi)));     % fine mesh size
ch = 1./(10*round(sqrt(2./fh)/10));         % coarse mesh size

% wavelength
high_wl = 2*pi./high_omega;
low_wl = 2*pi./low_omega;

% width of PML
high_wpml = 0.075*ones(size(high_omega));
low_wpml = 0.3*ones(size(high_omega));


%% Generate the domain sizes
sd = 1/2;
Rest = 0.35;           % estimate of the distance to the source point

high_r = NMLA_radius(high_omega(1),Rest);
md = sd + high_r + high_wpml;
md = ceil(md*10)/10;      % middle domain size 

% Rest = sqrt(2)*md;
low_r = NMLA_radius(low_omega(1),Rest);
ld = md + low_r + low_wpml;
ld = ceil(ld*10)/10;      % large domain size


%% Tests
tstart = tic;
for ti = 1: test_num
    omega = high_omega(ti);
    h = fh(ti);  h_c = ch(ti);
    fprintf(['-'*ones(1,80) '\n']);
    fprintf('\nCase %d: \nomega/(2*pi) = %d,   1/h = %d   1/h_c = %d,  NPW = %d ',...
        ti, round(omega/(2*pi)), 1/h, 1/h_c, NPW);
    
    
    %% Step 1: Solve the Hemholtz equation with the same source but with a relative low frequency sqrt(\omega) by Standard FEM, mesh size \omega*h = constant
    fprintf(['\n' '-'*ones(1,80) '\n']);
    fprintf('Step1: S-FEM, low frequency \n');
    tic;
    omega = low_omega(ti);              % low frequency
    a = ld(ti);                         % large domain 
    wpml = low_wpml(ti);                % width of PML
    sigmaMax = 25/wpml;                 % Maximun absorbtion
    [lnode,lelem] = squaremesh([-a,a,-a,a],h);
    
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
    
    % exact solution of smooth part
    uex = (1-cf).*ub;
    uex(rr <= epsilon) = 0;
    
    % errors
    du = u_std - uex;
    idx = find( ~( (x <= max(x)-wpml).*(x >= min(x)+wpml)...
        .*(y <= max(y)-wpml).*(y >= min(y)+wpml) ) ); % index on PML
    du(idx) = 0;  uex(idx) = 0;
    
    sfem_max_err(ti) = norm(du,inf);
    sfem_rel_max_err(ti) = norm(du,inf)/norm(uex,inf);
    sfem_l2_err(ti) = norm(du)*h;
    sfem_rel_l2_err(ti) = norm(du)/norm(uex);
    toc;
end


%% 
figure(2);
subplot(2,2,1);
show_convergence_rate(high_omega(1:test_num),sfem_max_err(1:test_num),'omega','Max ');
subplot(2,2,2);
show_convergence_rate(high_omega(1:test_num),sfem_rel_max_err(1:test_num),'omega','Rel max ');
subplot(2,2,3);
show_convergence_rate(high_omega(1:test_num),sfem_l2_err(1:test_num),'omega','L2 ');
subplot(2,2,4);
show_convergence_rate(high_omega(1:test_num),sfem_rel_l2_err(1:test_num),'omega','Rel L2 ');



%% print results
fprintf( ['\n' '-'*ones(1,80) '\n']);
fprintf( 'omega:                   ');
fprintf( '&  %.2e  ',high_omega );
fprintf( '\nomega/2pi:               ');
fprintf( '&  %.2e  ',high_omega/(2*pi) );
fprintf( '\n\nGrid size h:             ');
fprintf( '&  %.2e  ',fh);
fprintf( '\n1/h:                     ');
fprintf( '&  %.2e  ',1./fh);

fprintf( ['\n' '-'*ones(1,80) '\n']);
fprintf( '\n\nL2 error:                ');
fprintf( '&  %1.2d  ',sfem_l2_err);
fprintf( '\n\nRelative L2 error:       ');
fprintf( '&  %1.2d  ',sfem_rel_l2_err);

fprintf( ['\n' '-'*ones(1,80) '\n']);


