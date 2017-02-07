%% Convergence test for homogenenous case with exact ray information

% add path
clear;
addpath(genpath('../../ifem/'));
addpath('../Methods/');
addpath('../NMLA/');
addpath('../Helmholtz_data/');
addpath('../Plots_Prints/');

% set up
pde = Helmholtz_data_four_point_sources;    % pde data
xs = 20; ys = 20;       % source location
Nray = 4;
fquadorder = 9;    % the order of numerical quadrature
solver = 'DIR';      % the type of solver for the resulting linear system
plt = 0;                   % plot the solution or not
sec_opt = 0;           % not using second order correction of NMLA

test_num = 4;             % we test test_num examples
NPW = 4;                   % number of points per wavelength

% frequency
high_omega = [40 80 120 160 240 320 400 480 640]*pi;
low_omega = 2*sqrt(high_omega);

% error
low_l2_rayerr = 0*high_omega;      % L_2 ray error of low-freq waves
high_l2_rayerr = 0*high_omega;     % L_2 ray error of high-freq waves

rel_l2_err_numray = 0*high_omega;   % relative L_2 error of the numerical Ray-FEM solution
rel_l2_err_exray = 0*high_omega;   % relative L_2 error of the exact Ray-FEM solution



% wavelength
high_wl = 2*pi./high_omega;
low_wl = 2*pi./low_omega;

% mesh size
fh = 1./(NPW*round(high_omega/(2*pi)));      % fine mesh size
ch = 1./(20*round(low_omega/(4*pi)));        % coarse mesh size


%% Generate the domain sizes
sd = 1/2;
Rest = 20;          

high_r = NMLA_radius(high_omega,Rest);
md = sd + high_r;
md = ceil(md*10)/10;      % middle domain size 

low_r = NMLA_radius(low_omega,Rest);
ld = md + low_r;
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
    [lnode,lelem] = squaremesh([-a,a,-a,a],h);
    u_std = Standard_FEM_IBC(lnode,lelem,omega,pde,fquadorder,solver,plt,xs,ys);
    toc;
    
    
    %% Step 2: Use NMLA to find ray directions d_c with low frequency sqrt(\omega)
    fprintf([ '-'*ones(1,80) '\n']);
    fprintf('Step2: NMLA, low frequency \n');
    
    % compute numerical derivatives 
    [ux,uy] = num_derivative(u_std,h,2);
    
    a = md(ti);
    [mnode,melem] = squaremesh([-a,a,-a,a],h);
    
    [cnode,celem] = squaremesh([-a,a,-a,a],h_c);
    cN = size(cnode,1);    
    cnumray_angle = zeros(cN,Nray);
    
    % NMLA
    tic;
    for i = 1:cN
        x0 = cnode(i,1);  y0 = cnode(i,2);
        r0 = sqrt((x0-xs)^2 + (y0-ys)^2);
        c0 = pde.speed(cnode(i,:));
        [cnumray_angle(i,:)] = NMLA(x0,y0,c0,omega,Rest,lnode,lelem,u_std,ux,uy,[],1/8,Nray,'num',sec_opt,plt);
    end
    cnumray = exp(1i*cnumray_angle);
    numray1 = interpolation(cnode,celem,mnode,cnumray);
    toc;
    
    % compute the ray errors
    exray = pde.ray(mnode,xs,ys);
    rayerr1 = numray1 - exray;
    low_l2_rayerr(ti) = norm(rayerr1(:))*h/(norm(exray(:))*h);
    
    
    %% Step 3: Solve the original Helmholtz equation by Ray-based FEM with ray directions d_c
    fprintf(['-'*ones(1,80) '\n']);
    fprintf('Step3: Ray-FEM, high frequency \n');
   
    tic;
    omega = high_omega(ti);
    ray = numray1;
    uh1= Ray_FEM_IBC(mnode,melem,omega,pde,ray,fquadorder,plt,xs,ys);
    toc;
    
    
    
    %% Step 4: NMLA to find original ray directions d_o with wavenumber k
    fprintf(['-'*ones(1,80) '\n']);
    fprintf('Step4: NMLA, high frequency \n');
    
    % compute numerical derivatives
    [ux,uy] = num_derivative(uh1,h,2);
    
    a = sd;
    [node,elem] = squaremesh([-a,a,-a,a],h);
    
    [cnode,celem] = squaremesh([-a,a,-a,a],h_c);
    cN = size(cnode,1);
    cnumray_angle = zeros(cN,Nray);
    
    
    % NMLA
    tic;
    for i = 1:cN
        x0 = cnode(i,1);  y0 = cnode(i,2);
        r0 = sqrt((x0-xs)^2 + (y0-ys)^2);
        c0 = pde.speed(cnode(i,:));
        [cnumray_angle(i,:)] = NMLA(x0,y0,c0,omega,Rest,mnode,melem,uh1,ux,uy,[],1/8,Nray,'num',sec_opt,plt);
    end
    cnumray = exp(1i*cnumray_angle);
    numray2 = interpolation(cnode,celem,node,cnumray);
    toc;
    
    % compute the ray errors
    exray = pde.ray(node,xs,ys);
    rayerr2 = numray2 - exray;
    high_l2_rayerr(ti) = norm(rayerr2(:))*h/(norm(exray(:))*h);    
    
    
    %% Step 5: Solve the original Helmholtz equation by Ray-based FEM with ray directions d_o
    fprintf(['-'*ones(1,80) '\n']);
    fprintf('Step5: Ray-FEM, high frequency \n');
    tic;
    
    % Ray-FEM solution
    omega = high_omega(ti);
    ray = numray2;
    [uh2,~,~,~,rel_l2_err]= Ray_FEM_IBC(node,elem,omega,pde,ray,fquadorder,plt,xs,ys);
    rel_l2_err_numray(ti) = rel_l2_err;
    
    %% Ray-FEM with exact ray
    ray = pde.ray(node,xs,ys);
    [uh,~,~,~,rel_l2_err]= Ray_FEM_IBC(node,elem,omega,pde,ray,fquadorder,plt,xs,ys);
    rel_l2_err_exray(ti) = rel_l2_err;
    toc;
    
end

totaltime = toc(tstart);
fprintf('\n\nTotal running time: % d minutes \n', totaltime/60);



%% plots
figure(11);
subplot(1,2,1);
show_convergence_rate(high_omega(1:test_num),low_l2_rayerr(1:test_num),'omega','low rel l2');
subplot(1,2,2);
show_convergence_rate(high_omega(1:test_num),high_l2_rayerr(1:test_num),'omega','high rel l2');

figure(12);
subplot(1,2,1);
show_convergence_rate(high_omega(1:test_num),rel_l2_err_numray(1:test_num),'omega','NR-FEM rel l2');
subplot(1,2,2);
show_convergence_rate(high_omega(1:test_num),rel_l2_err_exray(1:test_num),'omega','ER-FEM rel l2');


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
fprintf( '\nLow rel l2 ray error:    ');
fprintf( '&  %1.2d  ',low_l2_rayerr);
fprintf( '\n\nHigh rel l2 ray error:   ');
fprintf( '&  %1.2d  ',high_l2_rayerr);

fprintf( '\n\nNR-FEM rel l2 error:    ');
fprintf( '&  %1.2d  ',rel_l2_err_numray);
fprintf( '\n\nER-FEM rel l2 error:     ');
fprintf( '&  %1.2d  ',rel_l2_err_exray);


fprintf( ['\n' '-'*ones(1,80) '\n']);


