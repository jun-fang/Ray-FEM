%% Convergence test for homogenenous case with numerical rays

clear;
addpath(genpath('../../ifem/'));
addpath('../Methods/');
addpath('../NMLA/');
addpath('../Cutoff_Functions/')
addpath('../Plots_Prints/');

%% Load source data
xs = 0; ys = 0;                     % point source location
epsilon = 50/(80*pi);                % cut-off parameter   
speed = @(x) ones(size(x,1),1);      % medium speed


%% Set up
plt = 0;                   % show solution or not
fquadorder = 3;            % numerical quadrature order
Nray = 1;                  % one ray direction at each grid node
sec_opt = 0;               % NMLA second order correction or not

NPW = 4;                   % number of points per wavelength
test_num = 5;              % we test test_num examples

% frequency
high_omega = [120 160 200 240 320 400 480 640]*pi;
low_omega = sqrt(high_omega); 

% error
low_max_rayerr = 0*high_omega;     % L_inf ray error of low-freq waves
low_l2_rayerr = 0*high_omega;      % L_2 ray error of low-freq waves
high_max_rayerr = 0*high_omega;    % L_inf ray error of high-freq waves
high_l2_rayerr = 0*high_omega;     % L_2 ray error of high-freq waves

max_err = 0*high_omega;            % L_inf error of the numerical solution
rel_max_err = 0*high_omega;        % relative L_inf error of the numerical solution
l2_err = 0*high_omega;             % L_2 error of the numerical solution
rel_l2_err = 0*high_omega;         % relative L_2 error of the numerical solution


% wavelength
high_wl = 2*pi./high_omega;
low_wl = 2*pi./low_omega;

% mesh size
fh = 1./(NPW*round(high_omega/(2*pi)));      % fine mesh size
ch = 1./(10*round(low_omega/(2*pi)));        % coarse mesh size
% ch = 1./(NPW*round(1./sqrt(fh)/10)*10);
% ch = fh.*ceil(ch./fh);

% width of PML
high_wpml = 0.15*ones(size(high_omega));
low_wpml = 0.35*ones(size(high_omega));
% high_wpml = 8*high_wl(1)*ones(size(high_omega)); %fh.*ceil(high_wl./fh);
% low_wpml = ch.*ceil(low_wl(1)./ch);


%% Generate the domain sizes
sd = 1/2;
Rest = 0.75;           % estimate of the distance to the source point

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
    A = assemble_Helmholtz_matrix_SFEM(lnode,lelem,omega,wpml,sigmaMax,speed,fquadorder);
    b = assemble_RHS_SFEM(lnode,lelem, @(x)nodal_basis(xs,ys,h,x),fquadorder);
    b = b/(h*h/2);
    [~,~,isBdNode] = findboundary(lelem);
    freeNode = find(~isBdNode);
    lN = size(lnode,1);        u_std = zeros(lN,1);
    u_std(freeNode) = A(freeNode,freeNode)\b(freeNode);
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
        c0 = speed(cnode(i,:));
        if r0 > (2*epsilon - 2*h)
            Rest = r0;
            [cnumray_angle(i)] = NMLA(x0,y0,c0,omega,Rest,lnode,lelem,u_std,ux,uy,[],1/5,Nray,'num',sec_opt,plt);
        else
            cnumray_angle(i) = ex_ray([x0,y0],xs,ys,0);
        end
    end
    cnumray = exp(1i*cnumray_angle);
    numray1 = interpolation(cnode,celem,mnode,cnumray);
    numray1 = numray1./abs(numray1);
    toc;
    
    % compute the ray errors
    exray = ex_ray(mnode,xs,ys,1);
    mr = sqrt((mnode(:,1)-xs).^2 + (mnode(:,2)-ys).^2);
    numray1 = numray1.*(mr > 2*epsilon) + exray.*(mr <= 2*epsilon);
    rayerr1 = numray1 - exray;
    low_max_rayerr(ti) = norm(rayerr1,inf);
    low_l2_rayerr(ti) = norm(rayerr1)*h/(norm(exray)*h);
    
    
    %% Step 3: Solve the original Helmholtz equation by Ray-based FEM with ray directions d_c
    fprintf(['-'*ones(1,80) '\n']);
    fprintf('Step3: Ray-FEM, high frequency \n');
   
    tic;
    omega = high_omega(ti);
    wpml = high_wpml(ti);                % width of PML
    sigmaMax = 25/wpml;                 % Maximun absorbtion
    ray = numray1;

    % smooth part
    option = 'homogeneous';
    A = assemble_Helmholtz_matrix_RayFEM(mnode,melem,omega,wpml,sigmaMax,speed,ray,fquadorder);
    b = assemble_RHS_RayFEM_with_ST(mnode,melem,xs,ys,omega,epsilon,wpml,sigmaMax,ray,speed,fquadorder,option);
    uh = RayFEM_direct_solver(mnode,melem,A,b,omega,ray,speed);
    
    % singularity part
    x = mnode(:,1); y = mnode(:,2);
    rr = sqrt((x-xs).^2 + (y-ys).^2);
    ub = 1i/4*besselh(0,1,omega*rr);
    cf = cutoff(epsilon,2*epsilon,mnode,xs,ys);
    
    % smooth + singularity
    uh1 = uh + ub.*cf;
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
        c0 = speed(cnode(i,:));
        if r0 > (2*epsilon - 2*h)
            Rest = r0;
            [cnumray_angle(i)] = NMLA(x0,y0,c0,omega,Rest,mnode,melem,uh1,ux,uy,[],1/5,Nray,'num',sec_opt,plt);
        else
            cnumray_angle(i) = ex_ray([x0,y0],xs,ys,0);
        end
    end
    cnumray = exp(1i*cnumray_angle);
    numray2 = interpolation(cnode,celem,node,cnumray);
    numray2 = numray2./abs(numray2);
    toc;
    
    % compute the ray errors
    exray = ex_ray(node,xs,ys,1);
    sr = sqrt((node(:,1)-xs).^2 + (node(:,2)-ys).^2);
    numray2 = numray2.*(sr > 2*epsilon) + exray.*(sr <= 2*epsilon);
    rayerr2 = numray2 - exray;
    high_max_rayerr(ti) = norm(rayerr2,inf);
    high_l2_rayerr(ti) = norm(rayerr2)*h/(norm(exray)*h);    
    
    
    %% Step 5: Solve the original Helmholtz equation by Ray-based FEM with ray directions d_o
    fprintf(['-'*ones(1,80) '\n']);
    fprintf('Step5: Ray-FEM, high frequency \n');
    tic;
    
    % Ray-FEM solution
    omega = high_omega(ti);
    wpml = high_wpml(ti);                % width of PML
    sigmaMax = 25/wpml;                 % Maximun absorbtion
    ray = numray2;
    
    option = 'homogeneous';
    A = assemble_Helmholtz_matrix_RayFEM(node,elem,omega,wpml,sigmaMax,speed,ray,fquadorder);
    b = assemble_RHS_RayFEM_with_ST(node,elem,xs,ys,omega,epsilon,wpml,sigmaMax,ray,speed,fquadorder,option);
    u = RayFEM_direct_solver(node,elem,A,b,omega,ray,speed);
    toc;
    
    % Excat solution 
    x = node(:,1); y = node(:,2);
    rr = sqrt((x-xs).^2 + (y-ys).^2);
    
    ub = 1i/4*besselh(0,1,omega*rr);
    cf = cutoff(epsilon,2*epsilon,node,xs,ys);
    uex = (1-cf).*ub;
    uex(rr <= epsilon) = 0;
   
    % Errors
    du = u - uex;
    idx = find( ~( (x <= max(x)-wpml).*(x >= min(x)+wpml)...
        .*(y <= max(y)-wpml).*(y >= min(y)+wpml) ) ); % index on PML
    du(idx) = 0;  uex(idx) = 0;
    
    max_err(ti) = norm(du,inf);
    rel_max_err(ti) = norm(du,inf)/norm(uex,inf);
    l2_err(ti) = norm(du)*h;
    rel_l2_err(ti) = norm(du)/norm(uex);
       
end

totaltime = toc(tstart);
fprintf('\n\nTotal running time: % d minutes \n', totaltime/60);



%% save output 
nameFile = strcat('resutls_2_HomNumRay_NPW_', num2str(NPW), '.mat');
save(nameFile, 'rel_l2_err', 'NPW', 'high_omega');



%% plots
figure(1);
subplot(2,2,1);
show_convergence_rate(high_omega(1:test_num),low_max_rayerr(1:test_num),'omega','low max');
subplot(2,2,2);
show_convergence_rate(high_omega(1:test_num),low_l2_rayerr(1:test_num),'omega','low l2');
subplot(2,2,3);
show_convergence_rate(high_omega(1:test_num),high_max_rayerr(1:test_num),'omega','high max');
subplot(2,2,4);
show_convergence_rate(high_omega(1:test_num),high_l2_rayerr(1:test_num),'omega','high l2');

figure(2);
subplot(2,2,1);
show_convergence_rate(high_omega(1:test_num),max_err(1:test_num),'omega','max err');
subplot(2,2,2);
show_convergence_rate(high_omega(1:test_num),l2_err(1:test_num),'omega','L2 err');
subplot(2,2,3);
show_convergence_rate(high_omega(1:test_num),rel_max_err(1:test_num),'omega','Rel max ');
subplot(2,2,4);
show_convergence_rate(high_omega(1:test_num),rel_l2_err(1:test_num),'omega','Rel L2 ');


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
fprintf( 'Low max ray error:       ');
fprintf( '&  %1.2d  ',low_max_rayerr);
fprintf( '\n\nLow rel l2 ray error:    ');
fprintf( '&  %1.2d  ',low_l2_rayerr);
fprintf( '\n\nHigh max ray error:      ');
fprintf( '&  %1.2d  ',high_max_rayerr);
fprintf( '\n\nHigh rel l2 ray error:   ');
fprintf( '&  %1.2d  ',high_l2_rayerr);
fprintf( '\n\nMax error:               ');
fprintf( '&  %1.2d  ',max_err);
fprintf( '\n\nRelative max error:      ');
fprintf( '&  %1.2d  ',rel_max_err);
fprintf( '\n\nL2 error:                ');
fprintf( '&  %1.2d  ',l2_err);
fprintf( '\n\nRelative L2 error:       ');
fprintf( '&  %1.2d  ',rel_l2_err);


fprintf( ['\n' '-'*ones(1,80) '\n']);


