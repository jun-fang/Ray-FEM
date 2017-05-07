%% h convergence test for homogenenous case with exact ray information
function [] = test6_cons_grad_velo(NPW, test_num)
% clear;
addpath(genpath('../../ifem/'));
addpath('../Methods/');
addpath('../Functions/')
addpath('../Plots_Prints/');


%% Set up
xs = 0; ys = 0;                     % source location
epsilon = 1/(2*pi);                 % cut-off parameter

v0 = 1; G0 = [0.1, -0.2];
speed = @(p) ( v0 + (G0(1)*p(:,1) + G0(2)*p(:,2)) );    % wave speed
speed_min = 0.85;

plt = 0;                   % show solution or not
fquadorder = 3;            % numerical quadrature order
Nray = 1;                  % one ray direction
sec_opt = 0;               % NMLA second order correction or not

high_omega = [120 160 240 320 400 500 600]*pi;
low_omega = 2*sqrt(high_omega);

% NPW = 4;
% test_num = 3;

% error
max_err = 0*zeros(1,test_num);      % L_inf error of the numerical solution
rel_max_err = 0*zeros(1,test_num);  % relative L_inf error of the numerical solution
l2_err = 0*zeros(1,test_num);       % L_2 error of the numerical solution
rel_l2_err = 0*zeros(1,test_num);   % relative L_2 error of the numerical solution
ref_l2 = 0*zeros(1,test_num);       % reference l2 norm

% wavelength
high_wl = 2*pi*speed_min./high_omega;
low_wl = 2*pi*speed_min./low_omega;

% mesh size
fh = 1./(10*round(NPW*high_omega/(2*pi*speed_min)/10));      % fine mesh size
ch = 1./(10*round(sqrt(4./fh)/10));                    % coarse mesh size

% width of PML
high_wpml = 0.07*ones(size(high_omega));
low_wpml = 0.19*ones(size(high_omega));


%% Generate the domain sizes
sd = 1/2;
Rest = 0.4654;            % estimate of the distance to the source point

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
    fprintf('\ncase %d: \nomega/(2*pi) = %d,   1/h = %d   1/h_c = %d,  NPW = %d \n',...
        ti, round(omega/(2*pi)), 1/h, 1/h_c, NPW);
    
    
    %% load Babich expansion
    switch round(omega/(pi))
        case 120
            load('Babich_CGV_30.mat');
        case 160
            load('Babich_CGV_40.mat');
        case 240
            load('Babich_CGV_60.mat');
        case 320
            load('Babich_CGV_80.mat');
        case 400
            load('Babich_CGV_100.mat');
        case 500
            load('Babich_CGV_125.mat');
        case 600
            load('Babich_CGV_150.mat');
    end
    
    a = sd;  Bx = -a: h : a;  By = -a: h : a;
    [BX0, BY0] = meshgrid(Bx0, By0);
    [BX, BY] = meshgrid(Bx, By);
    
    % refined phase and amplitude
    DD1 = interp2(BX0,BY0,D1,BX,BY,'spline'); % refined amplitude
    DD2 = interp2(BX0,BY0,D2,BX,BY,'spline'); % refined amplitude
    
    ttao = interp2(BX0,BY0,tao,BX,BY,'spline'); % refined phase
    mid = round(size(tao,1)/2);
    taox = tao2x ./ (2*tao);   taox(mid, mid) = 0;
    taoy = tao2y ./ (2*tao);   taoy(mid, mid) = 0;
    ttaox = interp2(BX0,BY0,taox,BX,BY,'spline'); % refined phase
    ttaoy = interp2(BX0,BY0,taoy,BX,BY,'spline'); % refined phase
    
    ray = atan2(ttaoy(:),ttaox(:));
    xx = BX(:)-xs;  yy = BY(:)-ys;
    rr = sqrt(xx.^2 + yy.^2);
    Bray = exp(1i*ray).*(rr>10*eps);
    Bray = reshape(Bray, size(BX));
    
    
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
    option ={ 'Babich', 'CGV', 'numerical_phase', high_omega(ti)};
    A = assemble_Helmholtz_matrix_SFEM(lnode,lelem,omega,wpml,sigmaMax,speed,fquadorder);
    b = assemble_RHS_SFEM_with_ST(lnode,lelem,xs,ys,omega,wpml,sigmaMax,epsilon,fquadorder,option);
    [~,~,isBdNode] = findboundary(lelem);
    freeNode = find(~isBdNode);
    lN = size(lnode,1);        u_std = zeros(lN,1);
    u_std(freeNode) = A(freeNode,freeNode)\b(freeNode);
    
    % singular part
    G1 = 1i*(sqrt(pi)/2)*besselh(0,1,omega*ttao);
    G1 = G1.*DD1;
    
    G2 = 1i*(sqrt(pi)/2)*exp(1i*pi)*(2*ttao/omega);
    G2 = G2.*besselh(1,1,omega*ttao);
    G2 = G2.*DD2;
    
    ln = round(sqrt(lN));
    ub = zeros(ln, ln);
    idx = round((ld-sd)/h)+1:round((ld+sd)/h)+1;
    ub(idx, idx) = G1 + G2;
    ub = ub(:);
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
            [cnumray_angle(i)] = NMLA(x0,y0,c0,omega,Rest,lnode,lelem,u_low,ux,uy,[],1/5,Nray,'num',sec_opt,plt);
        else
            cnumray_angle(i) =  ex_ray([x0,y0],xs,ys,0);
        end
    end
    cnumray = exp(1i*cnumray_angle);
    ray = interpolation(cnode,celem,mnode,cnumray);
    ray = ray./abs(ray);
    
    % Babich ray directions in the support of the cut-off dunction
    x = mnode(:,1); y = mnode(:,2);
    rr = sqrt((x-xs).^2 + (y-ys).^2);
    mN = length(rr); mn = round(sqrt(mN));
    mray = zeros(mn, mn);
    idx = round((md-sd)/h)+1:round((md+sd)/h)+1;
    mray(idx, idx) = Bray;
    mray = mray(:);
    ray(rr <= 2*epsilon) = mray(rr <= 2*epsilon);
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
    option ={ 'Babich', 'CGV', 'numerical_phase'};
    A = assemble_Helmholtz_matrix_RayFEM(mnode,melem,omega,wpml,sigmaMax,speed,ray,fquadorder);
    b = assemble_RHS_RayFEM_with_ST(mnode,melem,xs,ys,omega,epsilon,wpml,sigmaMax,ray,speed,fquadorder,option);
    uh = RayFEM_direct_solver(mnode,melem,A,b,omega,ray,speed);
    
    % singularity part
    G1 = 1i*(sqrt(pi)/2)*besselh(0,1,omega*ttao);
    G1 = G1.*DD1;
    
    G2 = 1i*(sqrt(pi)/2)*exp(1i*pi)*(2*ttao/omega);
    G2 = G2.*besselh(1,1,omega*ttao);
    G2 = G2.*DD2;
    
    ub = zeros(mn, mn);
    idx = round((md-sd)/h)+1:round((md+sd)/h)+1;
    ub(idx, idx) = G1 + G2;
    ub = ub(:);
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
            [cnumray_angle(i)] = NMLA(x0,y0,c0,omega,Rest,mnode,melem,uh1,ux,uy,[],1/5,Nray,'num',sec_opt,plt);
        else
            cnumray_angle(i) =  ex_ray([x0,y0],xs,ys,0);
        end
    end
    cnumray = exp(1i*cnumray_angle);
    ray = interpolation(cnode,celem,node,cnumray);
    ray = ray./abs(ray);
    
    % Babich ray directions in the support of the cut-off dunction
    x = node(:,1); y = node(:,2);
    rr = sqrt((x-xs).^2 + (y-ys).^2);
    sray = Bray(:);
    ray(rr <= 2*epsilon) = sray(rr <= 2*epsilon);
    toc;
    
    
    
    %% Step 5: Solve the original Helmholtz equation by Ray-based FEM with ray directions d_o
    fprintf([ '-'*ones(1,80) '\n']);
    fprintf('Step5: Ray-FEM, high frequency \n');
    tic;
    
    % Assembling
    option ={ 'Babich', 'CGV', 'numerical_phase'};
    A = assemble_Helmholtz_matrix_RayFEM(node,elem,omega,wpml,sigmaMax,speed,ray,fquadorder);
    b = assemble_RHS_RayFEM_with_ST(node,elem,xs,ys,omega,epsilon,wpml,sigmaMax,ray,speed,fquadorder,option);
    uh = RayFEM_direct_solver(node,elem,A,b,omega,ray,speed);
    toc;
    
    
    %% Compute errors
    fprintf([ '-'*ones(1,80) '\n']);
    fprintf('Compute errors\n');
    tic;
    
    % reference solution
    ub = G1 + G2; ub = ub(:);
    cf = cutoff(epsilon,2*epsilon,node,xs,ys);
    ur = (1-cf).*ub;
    ur(rr <= epsilon) = 0;
    
    % Errors
    du = uh - ur;
    idx = find( ~( (x<=max(x)-wpml).*(x>= min(x)+wpml)...
        .*(y<= max(y)-wpml).*(y>= min(y)+wpml) ) ); % index on PML
    du(idx) = 0;  ur(idx) = 0;
    
    max_err(ti) = norm(du,inf);
    rel_max_err(ti) = norm(du,inf)/norm(ur,inf);
    l2_err(ti) = norm(du)*h;
    ref_l2(ti) = norm(ur)*h;
    rel_l2_err(ti) = l2_err(ti)/ref_l2(ti)
    
    toc;
    
    
end

totaltime = toc(tstart);
fprintf('\n\nTotal running time: % d minutes \n', totaltime/60);

nameFile = strcat('resutls_6_CGV_NPW_', num2str(NPW), '.mat');
save(nameFile, 'rel_l2_err', 'NPW', 'high_omega', 'test_num');

figure(62);
show_convergence_rate(high_omega(1:test_num),rel_l2_err(1:test_num),'h','Rel L2 err');



