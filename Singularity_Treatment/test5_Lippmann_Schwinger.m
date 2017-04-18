%% Test for Lippmann-Schwinger equation

clear;

addpath(genpath('../../ifem/'));
addpath('../Methods/');
addpath('../NMLA/');
addpath('../Cutoff_Functions/')
addpath('../Plots_Prints/');
% addpath('/home/jun/Documents/MATLAB/Solutions_Lippmann_Schwinger');


%% Load source/wavespeed data
xs = -0.2; ys = -0.2;                     % point source location

sigma = 0.15;
xHet = 0.2;   yHet = 0.2;

nu = @(x,y) 0.2*exp( -1/(2*sigma^2)*((x-xHet).^2 + (y-yHet).^2) )...
    .*Lippmann_Schwinger_window(sqrt((x-xHet).^2 + (y-yHet).^2), 0.22,0.16  );

speed = @(p) 1./sqrt(1 + nu( p(:,1), p(:,2) ));    % wave speed
speed_min = 1/sqrt(1.2);


%% Set up
plt = 0;                   % show solution or not
fquadorder = 3;            % numerical quadrature order
Nray = 1;                  % one ray direction
sec_opt = 0;               % NMLA second order correction or not
epsilon = 1/(2*pi);        % cut-off parameter


NPW = 4;                   % number of points per wavelength
test_num = 1;              % we test test_num examples

% frequency
high_omega = [120 160 240 320 480 640 800 960]*pi;
high_omega = high_omega(3:end);
low_omega = 2*sqrt(high_omega);

% error
max_err = 0*high_omega;            % L_inf error of the numerical solution
rel_max_err = 0*high_omega;        % relative L_inf error of the numerical solution
l2_err = 0*high_omega;             % L_2 error of the numerical solution
rel_l2_err = 0*high_omega;         % relative L_2 error of the numerical solution
ref_l2 = 0*high_omega;             % reference l2 norm

rays = cell(test_num,1);

% wavelength
high_wl = 2*pi*speed_min./high_omega;
low_wl = 2*pi*speed_min./low_omega;

% mesh size
fh = 1./(10*round(NPW*high_omega/(2*pi*speed_min)/10));      % fine mesh size
ch = 1./(10*round(sqrt(2./fh)/10));                    % coarse mesh size

% width of PML
high_wpml = 0.065*ones(size(high_omega));
low_wpml = 0.27*ones(size(high_omega));


%% Generate the domain sizes
sd = 1/2;
% Rest = 0.35;
% Rest = 0.5618;
Rest = sqrt((sd-xs)^2 + (sd-ys)^2);           % estimate of the distance to the source point

high_r = NMLA_radius(high_omega(1),Rest);
md = sd + high_r + high_wpml;
md = ceil(md*10)/10;      % middle domain size

Rest = sqrt((md(1)-xs)^2 + (md(1)-ys)^2);           % estimate of the distance to the source point
low_r = NMLA_radius(low_omega(1),Rest);
ld = md + low_r + low_wpml;
ld = ceil(ld*10)/10;      % large domain size


%% Tests
tstart = tic;
for ti = 1: test_num
    omega = high_omega(ti);
    h = fh(ti);  h_c = ch(ti);
    fprintf(['-'*ones(1,80) '\n']);
    fprintf('\ncase %d: \nomega/(2*pi) = %d,   1/h = %d   1/h_c = %d,  NPW = %d \n',...
        ti, round(omega/(2*pi)), 1/h, 1/h_c, NPW);
    
    
    %% Step 1: Solve the Hemholtz equation with the same source but with a relative low frequency sqrt(\omega) by Standard FEM, mesh size \omega*h = constant
    fprintf(['-'*ones(1,80) '\n']);
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
    
    % low frequency solution: smooth + singularity
    u_low = u_std + ub.*cf;
    toc;
    
    
    %% Step 2: Use NMLA to find ray directions d_c with low frequency sqrt(\omega)
    fprintf(['-'*ones(1,80) '\n']);
    fprintf('Step2: NMLA, low frequency \n');
    tic;
    
    % compute numerical derivatives
    [ux,uy] = num_derivative(u_low,h,2);
    
    a = md(ti);
    [mnode,melem] = squaremesh([-a,a,-a,a],h);
    [cnode,celem] = squaremesh([-a,a,-a,a],h_c);
    cN = size(cnode,1);
    cnumray_angle = zeros(cN,Nray);
    
    % NMLA
    for i = 1:cN
        x0 = cnode(i,1);  y0 = cnode(i,2);
        r0 = sqrt((x0-xs)^2 + (y0-ys)^2);
        c0 = speed(cnode(i,:));
        if r0 > (2*epsilon - 3*h_c)
            Rest = r0;
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
%     ray((x < xs) & (y < ys)) = exray((x < xs) & (y < ys));
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
    [ux,uy] = num_derivative(uh1,h,2);
    
    a = sd;
    [node,elem] = squaremesh([-a,a,-a,a],h);
    [cnode,celem] = squaremesh([-a,a,-a,a],h_c);
    cN = size(cnode,1);
    cnumray_angle = zeros(cN,Nray);
     
    % NMLA
    for i = 1:cN
        x0 = cnode(i,1);  y0 = cnode(i,2);
        r0 = sqrt((x0-xs)^2 + (y0-ys)^2);
        c0 = speed(cnode(i,:));
        if r0 > (2*epsilon - 3*h_c)
            Rest = r0;
            [cnumray_angle(i)] = NMLA(x0,y0,c0,omega,Rest,mnode,melem,uh1,ux,uy,[],1/5,Nray,'num',sec_opt,plt);
        else
            cnumray_angle(i) =  ex_ray([x0,y0],xs,ys,0);
        end
    end
    cnumray = exp(1i*cnumray_angle);
    ray = interpolation(cnode,celem,node,cnumray);
    ray = ray./abs(ray);
    
    % analytical ray directions in the support of the cut-off dunction
    exray = ex_ray(node,xs,ys,1);
    x = node(:,1); y = node(:,2);
    rr = sqrt((x-xs).^2 + (y-ys).^2);
    ray(rr <= 2*epsilon) = exray(rr <= 2*epsilon);
%     ray((x < xs) & (y < ys)) = exray((x < xs) & (y < ys));
%     figure(2); ray_field(ray,node,10,1/10);
    toc;
    
    
    %% Step 5: Solve the original Helmholtz equation by Ray-based FEM with ray directions d_o
    fprintf([ '-'*ones(1,80) '\n']);
    fprintf('Step5: Ray-FEM, high frequency \n');
    tic;
    
%     ray = exray;
    % Assembling
    omega = high_omega(ti);
    wpml = high_wpml(ti);                % width of PML
    sigmaMax = 25/wpml;                 % Maximun absorbtion
    
    option = 'homogeneous';
    A = assemble_Helmholtz_matrix_RayFEM(node,elem,omega,wpml,sigmaMax,speed,ray,fquadorder);
    b = assemble_RHS_RayFEM_with_ST(node,elem,xs,ys,omega,epsilon,wpml,sigmaMax,ray,speed,fquadorder,option);
    uh = RayFEM_direct_solver(node,elem,A,b,omega,ray,speed);
%     [~,v] = RayFEM_direct_solver(node,elem,A,b,omega,ray,speed);
    toc;
    
    clear lnode lelem mnode melem cnode celem cnumray cnumray_angle;
    clear A b x y rr exray cf u_std ub uh1 ux uy freeNode isBdNode;
    
    
    %% Compute errors
    fprintf([ '-'*ones(1,80) '\n']);
    fprintf('Compute errors\n');
    tic;
    switch round(omega/(pi))
        case 120
            load('/Solutions_Lippmann_Schwinger/point_source_k_120_jun_new_version.mat')
        case 160
            load('/Solutions_Lippmann_Schwinger/point_source_k_160_jun_new_version.mat')
        case 240
            load('/Solutions_Lippmann_Schwinger/point_source_k_240_jun_new_version.mat')
        case 320
            load('/Solutions_Lippmann_Schwinger/point_source_k_320_jun_new_version.mat')
        case 480
            load('/Solutions_Lippmann_Schwinger/point_source_k_480_jun_new_version.mat')
        case 640
            load('/Solutions_Lippmann_Schwinger/point_source_k_640_jun_new_version.mat')
        case 960
            load('/Solutions_Lippmann_Schwinger/point_source_k_960_jun_new_version.mat')
    end

    rh = 1/5000;
    [rnode,relem] = squaremesh([-a,a,-a,a],rh); rn = round(sqrt(size(rnode,1)));
    ur = interpolation(rnode,relem,node,u);
%     uh = RayFEM_solution(node,elem,omega,speed,v,ray,rnode);
    clear rnode relem;
    
    % Reference solution 
    rnode = node;
    x = rnode(:,1); y = rnode(:,2);
    rr = sqrt((x-xs).^2 + (y-ys).^2);
    
    cf = cutoff(epsilon,2*epsilon,rnode,xs,ys);
    ur = (1-cf).*ur;
    ur(rr<=epsilon) = 0;
   
    % Errors
    du = uh - ur;
    idx = find( ~( (x<=max(x)-wpml).*(x>= min(x)+wpml)...
        .*(y<= max(y)-wpml).*(y>= min(y)+wpml) ) ); % index on PML
%     idx = find( ~( (x<=0.4).*(x>= -0.2)...
%         .*(y<= 0.4).*(y>= -0.2) ) ); % index on PML
    du(idx) = 0;  ur(idx) = 0;
    
    max_err(ti) = norm(du,inf);
    rel_max_err(ti) = norm(du,inf)/norm(ur,inf);
    l2_err(ti) = norm(du)*h;
    ref_l2(ti) = norm(ur)*h;
    rel_l2_err(ti) = l2_err(ti)/ref_l2(ti);
        
    toc;
    
%     rays{ti} = ray;
%     clear node elem rnode x y rr cf idx u;
    
    
end

totaltime = toc(tstart);
fprintf('\n\nTotal running time: % d minutes \n', totaltime/60);


%% saving output 
nameFile = strcat('resutls_5_LipSch_NPW_', num2str(NPW), '.mat');
save(nameFile, 'ref_l2', 'l2_err' , 'rel_l2_err', 'high_omega', 'test_num');


%% plots


figure(8);
subplot(1,2,1);
showsolution(node,elem,real(du(:))); colorbar;
title('LipSch Ray-FEM error')
subplot(1,2,2);
showsolution(node,elem,real(du(:)),2); colorbar;
title('LipSch Ray-FEM error')

% subplot(2,2,3);
% showsolution(node,elem,real(uh(:))); colorbar;
% title('Ray-FEM solution error')
% subplot(2,2,4);
% showsolution(node,elem,real(uh(:)),2); colorbar;
% title('Ray-FEM solution error')


% figure(3);
% % reference L2 norm
% subplot(1,3,1);
% show_convergence_rate(high_omega(1:test_num),ref_l2(1:test_num),'omega','Ref L2');
% % Ray-FEM solution L2 error
% subplot(1,3,2);
% show_convergence_rate(high_omega(1:test_num),l2_err(1:test_num),'omega','L2 err');
% % Ray-FEM solution relative L2 error
% subplot(1,3,3);
% show_convergence_rate(high_omega(1:test_num),rel_l2_err(1:test_num),'omega','Rel L2 err');



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


% --------------------------------------------------------------------------------
% Max error:               &  5.35e-03  &  2.81e-03  &  3.09e-03  &  2.31e-03  &  1.83e-03  &  00  
% 
% Relative max error:      &  2.86e-01  &  1.61e-01  &  2.24e-01  &  1.99e-01  &  1.98e-01  &  00  
% 
% L2 error:                &  1.96e-03  &  1.09e-03  &  4.69e-04  &  2.53e-04  &  1.13e-04  &  00  
% 
% Relative L2 error:       &  1.31e-02  &  1.12e-02  &  8.83e-03  &  7.33e-03  &  6.02e-03  &  00  
% --------------------------------------------------------------------------------


%     % Plots
%     sh = 1/500;
%     [snode,selem] = squaremesh([-a,a,-a,a],sh); 
%     sr = sqrt((snode(:,1)-xs).^2 + (snode(:,2)-ys).^2);
%     sn = round(sqrt(size(snode,1))); idx = 1:sn; idx = 10*(idx-1) + 1;
%     su1 = reshape(du,rn,rn); su1 = su1(idx,idx);
%     su2 = reshape(ur,rn,rn); su2 = su2(idx,idx);
%     figure(8); 
%     subplot(2,2,1);
%     showsolution(snode,selem,real(su1(:))); colorbar;
%     title('Ray-FEM solution error')
%     subplot(2,2,2);
%     showsolution(snode,selem,real(su1(:)),2); colorbar;
%     title('Ray-FEM solution error')
%     subplot(2,2,3);
%     showsolution(snode,selem,real(su2(:))); colorbar;
%     title('Ray-FEM solution error')
%     subplot(2,2,4);
%     showsolution(snode,selem,real(su2(:)),2); colorbar;
%     title('Ray-FEM solution error')
    
%     figure(9);
%     subplot(1,2,1);
%     showsolution(snode,selem,speed(snode).*(sr>=2*epsilon)); colorbar;
%     title('Wave speed')
%     subplot(1,2,2);
%     showsolution(snode,selem,speed(snode).*(sr>=2*epsilon),2); colorbar;
%     title('Wave speed')