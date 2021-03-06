%% h convergence test for homogenenous case with numerical ray information
%    O(h^2)


clear;
addpath(genpath('../../../ifem/'));
addpath('../../Methods/');
addpath('../../NMLA/');
addpath('../../Functions/')
addpath('../../Plots/');


%% Load source/wavespeed data
xs = 0; ys = 0.35;               % point source location
speed = @(p) ones(size(p(:,1)));    % wave speed
speed_min = 1;


%% Set up
plt = 0;                   % show solution or not
fquadorder = 3;            % numerical quadrature order
Nray = 1;                  % one ray direction
sec_opt = 0;               % NMLA second order correction or not
epsilon = 1/(2.2*pi);        % cut-off parameter

high_omega = 100*pi;
low_omega = 4*sqrt(high_omega);

% width of PML
high_wpml = 0.05;
low_wpml = 0.1;

Rest = 1.5;
high_r = NMLA_radius(high_omega/3,Rest);
low_r = NMLA_radius(low_omega/3,Rest);

low_d = round((low_r + low_wpml)/0.05)*0.05;
high_d = round((high_r + high_wpml)/0.05)*0.05;


% wavelength
high_wl = 2*pi*speed_min./high_omega;
low_wl = 2*pi*speed_min./low_omega;

wl = high_wl;

% mesh size
fh = 1/800; %1/(10*round(16*high_omega/(2*pi*speed_min)/10));  % fine mesh size
ch = 1/80;                                             % coarse mesh size

% sd = 1/2; md = 0.6; ld = 0.9;

% domain
sdx = 1.5; sdy = 0.5;
mdx = sdx + high_d; mdy = sdy + high_d;
ldx = mdx + low_d; ldy = mdy + low_d;


large_domain = [-ldx, ldx, -ldy, ldy];
middle_domain = [-mdx, mdx, -mdy, mdy];
small_domain = [-sdx, sdx, -sdy, sdy];

% Rest = 0.5618; ti = 1;

Rest = 1; ti = 1;


h = fh; h_c = ch;
tstart = tic;

%% Step 1: Solve the Hemholtz equation with the same source but with a relative low frequency sqrt(\omega) by Standard FEM, mesh size \omega*h = constant
fprintf(['-'*ones(1,80) '\n']);
fprintf('Step1: S-FEM, low frequency \n');
tic;
omega = low_omega(ti);              % low frequency
% a = ld(ti);                         % large domain
wpml = low_wpml(ti);                % width of PML
sigmaMax = 25/wpml;                 % Maximun absorbtion
% [lnode,lelem] = squaremesh([-a,a,-a,a],h);
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


%% Step 2: Use NMLA to find ray directions d_c with low frequency sqrt(\omega)
fprintf(['-'*ones(1,80) '\n']);
fprintf('Step2: NMLA, low frequency \n');
tic;

m = round( (large_domain(2) - large_domain(1)) /h ) + 1;
n = round( (large_domain(4) - large_domain(3)) /h ) + 1;
[ux,uy] = num_derivative(u_low,h,2,m,n);

% % compute numerical derivatives
% [ux,uy] = num_derivative(u_low,h,2);

% a = md(ti);
% [mnode,melem] = squaremesh([-a,a,-a,a],h);
% [cnode,celem] = squaremesh([-a,a,-a,a],h_c);

[mnode,melem] = squaremesh(middle_domain,h);
[cnode,celem] = squaremesh(middle_domain,h_c);


cN = size(cnode,1);
cnumray_angle = zeros(cN,Nray);

% NMLA
for i = 1:cN
    x0 = cnode(i,1);  y0 = cnode(i,2);
    r0 = sqrt((x0-xs)^2 + (y0-ys)^2);
    c0 = speed(cnode(i,:));
    if r0 > (2*epsilon - 3*h_c)
        Rest = min(r0, 1);
        [cnumray_angle(i)] = NMLA(x0,y0,c0,omega,Rest,lnode,lelem,u_low,ux,uy,[],1/5,Nray,'num',sec_opt,plt);
    else
        cnumray_angle(i) =  ex_ray([x0,y0],xs,ys,0);
    end
end
cnumray = exp(1i*cnumray_angle);
ray = interpolation(cnode,celem,mnode,cnumray);
ray = ray./abs(ray);

% analytical ray directions in the support of the cut-off dunction
exray = ex_ray(mnode,xs,ys,1);
x = mnode(:,1); y = mnode(:,2);
rr = sqrt((x-xs).^2 + (y-ys).^2);
ray(rr <= 2*epsilon) = exray(rr <= 2*epsilon);
%     figure(1); ray_field(ray,mnode,10,1/10);
toc;


%% Step 3: Solve the original Helmholtz equation by Ray-based FEM with ray directions d_c
fprintf(['-'*ones(1,80) '\n']);
fprintf('Step3: Ray-FEM, high frequency \n');
tic;
omega = high_omega(ti);
wpml = high_wpml(ti);                % width of PML
sigmaMax = 25/wpml;                 % Maximun absorbtion

% smooth part
option = 'homogeneous';
A = assemble_Helmholtz_matrix_RayFEM(mnode,melem,omega,wpml,sigmaMax,speed,ray,fquadorder);
b = assemble_RHS_RayFEM_with_ST(mnode,melem,xs,ys,omega,epsilon,wpml,sigmaMax,ray,speed,fquadorder,option);
uh = RayFEM_direct_solver(mnode,melem,A,b,omega,ray,speed);

% singularity part
ub = 1i/4*besselh(0,1,omega*rr);
cf = cutoff(epsilon,2*epsilon,mnode,xs,ys);

% smooth + singularity
uh1 = uh + ub.*cf;
toc;



%% Step 4: NMLA to find original ray directions d_o with wavenumber k
fprintf(['-'*ones(1,80) '\n']);
fprintf('Step4: NMLA, high frequency \n');
tic;

% compute numerical derivatives
m = round( (middle_domain(2) - middle_domain(1)) /h ) + 1;
n = round( (middle_domain(4) - middle_domain(3)) /h ) + 1;
[ux,uy] = num_derivative(uh1,h,2,m,n);


% % compute numerical derivatives
% [ux,uy] = num_derivative(uh1,h,2);

% a = sd;
% [node,elem] = squaremesh([-a,a,-a,a],h);
% [cnode,celem] = squaremesh([-a,a,-a,a],h_c);

[node,elem] = squaremesh(small_domain,h);
[cnode,celem] = squaremesh(small_domain,h_c);

cN = size(cnode,1);
cnumray_angle = zeros(cN,Nray);

% NMLA
for i = 1:cN
    x0 = cnode(i,1);  y0 = cnode(i,2);
    r0 = sqrt((x0-xs)^2 + (y0-ys)^2);
    c0 = speed(cnode(i,:));
    if r0 > (2*epsilon - 3*h_c)
        Rest = min(r0, 1);
        [cnumray_angle(i)] = NMLA(x0,y0,c0,omega,Rest,mnode,melem,uh1,ux,uy,[],1/5,Nray,'num',sec_opt,plt);
    else
        cnumray_angle(i) =  ex_ray([x0,y0],xs,ys,0);
    end
end
cnumray = exp(1i*cnumray_angle);
toc;


% %% h refinement
% h_ray = 1./(10*round(32*high_omega/(2*pi*speed_min)/10));
% % [node_ray, elem_ray] = squaremesh([-a,a,-a,a],h_ray);
% [node_ray, elem_ray] = squaremesh(small_domain,h_ray);
% rray = interpolation(cnode,celem,node_ray,cnumray);

rray = interpolation(cnode,celem,node,cnumray);
rray = rray./abs(rray);

exray = ex_ray(node,xs,ys,1);
x = node(:,1); y = node(:,2);
rr = sqrt((x-xs).^2 + (y-ys).^2);
rray(rr <= 2*epsilon) = exray(rr <= 2*epsilon);

% Assembling
omega = high_omega;
wpml = high_wpml;                % width of PML
sigmaMax = 25/wpml;                 % Maximun absorbtion
ray = rray;

option = 'homogeneous';
A = assemble_Helmholtz_matrix_RayFEM(node,elem,omega,wpml,sigmaMax,speed,ray,fquadorder);
b = assemble_RHS_RayFEM_with_ST(node,elem,xs,ys,omega,epsilon,wpml,sigmaMax,ray,speed,fquadorder,option);
ur = RayFEM_direct_solver(node,elem,A,b,omega,ray,speed);
toc;

% singularity part
x = node(:,1); y = node(:,2);
rr = sqrt((x-xs).^2 + (y-ys).^2);
ub = 1i/4*besselh(0,1,omega*rr);
cf = cutoff(epsilon,2*epsilon,node,xs,ys);

% smooth + singularity
ur = ur + ub.*cf;

ur(rr<2*wl) = 0;
unorm = norm(ur)*h;



rh = h; rN = size(node,1);


rm = round( (small_domain(2) - small_domain(1)) /h ) + 1;
rn = round( (small_domain(4) - small_domain(3)) /h ) + 1;
ur = reshape(ur, rn, rm);
rray = reshape(rray, rn, rm);


NPW = [1,2,4, 8];
fh = 1./(10*round(NPW*high_omega/(2*pi*speed_min)/10));

test_num = 4;

% error
max_err = 0*zeros(1,test_num);      % L_inf error of the numerical solution
rel_max_err = 0*zeros(1,test_num);  % relative L_inf error of the numerical solution
l2_err = 0*zeros(1,test_num);       % L_2 error of the numerical solution
rel_l2_err = 0*zeros(1,test_num);   % relative L_2 error of the numerical solution
ref_l2 = 0*zeros(1,test_num);       % reference l2 norm



ti = 1;
%% tests
for ti = 1: test_num
    h = fh(ti);
    [node,elem] = squaremesh(small_domain,h);
    N = size(node,1);
    nh = round(h/rh);
    im = 1:nh:rm; in = 1:nh:rn;
    ray = rray(in, im);
    ray = ray(:);
    
    
   
    %% Step 5: Solve the original Helmholtz equation by Ray-based FEM with ray directions d_o
    fprintf([ '-'*ones(1,80) '\n']);
    fprintf('Step5: Ray-FEM, high frequency \n');
    tic;
    
    % Assembling
    omega = high_omega;
    wpml = high_wpml;                % width of PML
    sigmaMax = 25/wpml;                 % Maximun absorbtion
    
    option = 'homogeneous';
    A = assemble_Helmholtz_matrix_RayFEM(node,elem,omega,wpml,sigmaMax,speed,ray,fquadorder);
    b = assemble_RHS_RayFEM_with_ST(node,elem,xs,ys,omega,epsilon,wpml,sigmaMax,ray,speed,fquadorder,option);
    uh = RayFEM_direct_solver(node,elem,A,b,omega,ray,speed);
    toc;
    
    % singularity part
    x = node(:,1); y = node(:,2);
    rr = sqrt((x-xs).^2 + (y-ys).^2);
    ub = 1i/4*besselh(0,1,omega*rr);
    cf = cutoff(epsilon,2*epsilon,node,xs,ys);
    
    % smooth + singularity
    uh = uh + ub.*cf;

    
    %% Compute errors
    fprintf([ '-'*ones(1,80) '\n']);
    fprintf('Compute errors\n');
    tic;
    
%     rh = 1/5000;
%     [rnode,~] = squaremesh(small_domain,rh); rn = round(sqrt(size(rnode,1)));
% 
% %     [rnode,~] = squaremesh([-a,a,-a,a],rh); rn = round(sqrt(size(rnode,1)));
%     uh = RayFEM_solution(node,elem,omega,speed,v,ray,rnode);
%     
%     % Reference solution
%     x = rnode(:,1); y = rnode(:,2);
%     rr = sqrt((x-xs).^2 + (y-ys).^2);
%     
%     u = 1i/4*besselh(0,1,omega*rr);
%     cf = cutoff(epsilon,2*epsilon,rnode,xs,ys);
%     ur = (1-cf).*u;
%     ur(rr<=epsilon) = 0;
    
    % Errors
    urr = ur(in,im); urr = urr(:);
    uh(rr<2*wl) = 0;
    du = uh - urr;
    
    idx = find( ~( (x<=max(x)-wpml).*(x>= min(x)+wpml)...
        .*(y<= max(y)-wpml).*(y>= min(y)+wpml) ) ); % index on PML
    du(idx) = 0;  urr(idx) = 0;
    
    max_err(ti) = norm(du,inf);
    rel_max_err(ti) = norm(du,inf)/norm(urr,inf);
    l2_err(ti) = norm(du)*h;
    ref_l2(ti) = norm(ur(:))*rh;
    rel_l2_err(ti) = l2_err(ti)/ref_l2(ti)
    
    toc;

    
end

totaltime = toc(tstart);
fprintf('\n\nTotal running time: % d minutes \n', totaltime/60);

figure(22);
show_convergence_rate(fh(1:test_num),rel_l2_err(1:test_num),'h','Rel L2 err');


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

fprintf( '\n\nMax error:               ');
fprintf( '&  %1.2d  ',max_err);
fprintf( '\n\nRelative max error:      ');
fprintf( '&  %1.2d  ',rel_max_err);

fprintf( ['\n' '-'*ones(1,80) '\n']);

fprintf( '\n\nReference L2 norm:       ');
fprintf( '&  %1.2d  ',ref_l2);
fprintf( '\n\nL2 error:                ');
fprintf( '&  %1.2d  ',l2_err);
fprintf( '\n\nRelative L2 error:       ');
fprintf( '&  %1.2d  ',rel_l2_err);

fprintf( ['\n' '-'*ones(1,80) '\n']);

