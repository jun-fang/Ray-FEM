%% Phase error test for point source inside homogeneous medium domain 

% add path
clear;
addpath(genpath('../../ifem/'));
addpath('../Methods/');
addpath('../NMLA/');
addpath('../Helmholtz_data/');
addpath('../Plots_Prints/');

% set up
xs = -0.4; ys = -0.4;       % source location
speed = @(p) ones(size(p(:,1)));
Nray = 1;
fquadorder = 3;    % the order of numerical quadrature
solver = 'DIR';      % the type of solver for the resulting linear system
plt = 0;                   % plot the solution or not
sec_opt = 0;           % not using second order correction of NMLA

test_num = 1;             % we test test_num examples
NPWs = [6 8 10];                   % number of points per wavelength

% frequency
high_omega = 250*pi;
low_omega = sqrt(high_omega);


% wavelength
high_wl = 2*pi./high_omega;
low_wl = 2*pi./low_omega;

% source neighborhood size
r = 16*high_wl;

% width of PML
high_wpml = 4*high_wl(1)*ones(size(high_omega)); 
ch = 1./(20*round(low_omega/(4*pi)));        % coarse mesh size
low_wpml = ch.*ceil(low_wl(1)./ch);


%% Generate the domain sizes
sd = 1/2;
Rest = 1;           % estimate of the distance to the source point

high_r = NMLA_radius(high_omega,Rest);
md = sd + high_r + high_wpml;
md = ceil(md*10)/10;      % middle domain size

low_r = NMLA_radius(low_omega,Rest);
ld = md + low_r + low_wpml;
ld = ceil(ld*10)/10;      % large domain size


%% Tests
tstart = tic;
for ti = 1: test_num
    omega = high_omega;
    NPW = NPWs(ti);
    fh = 1./(NPW*round(high_omega/(2*pi)));      % fine mesh size
    ch = 1./(20*round(low_omega/(4*pi)));        % coarse mesh size
    
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
    cnumray = zeros(cN,Nray);
    
    % NMLA
    tic;
    for i = 1:cN
        x0 = cnode(i,1);  y0 = cnode(i,2);
        d0 = sqrt((x0-xs)^2 + (y0-ys)^2);
        c0 = speed(cnode(i,:));
        if d0 > r
            Rest = d0;
            [cnumray_angle] = NMLA(x0,y0,c0,omega,Rest,lnode,lelem,u_std,ux,uy,[],1/5,Nray,'num',sec_opt,plt);
            cnumray(i) = exp(1i*cnumray_angle);
        else
            cnumray(i) = ex_ray([x0,y0],xs,ys,1);
        end
    end
    
    numray1 = interpolation(cnode,celem,mnode,cnumray);
    toc;
    
    % compute the ray errors
    exray = ex_ray(mnode,xs,ys,1);
    xx = mnode(:,1) - xs;   yy = mnode(:,2) - ys;
    rr = sqrt(xx.*xx + yy.*yy);
    numray1(rr<=r) = exray(rr<=r);
    rayerr1 = numray1 - exray;
    
    
    %% Step 3: Solve the original Helmholtz equation by Ray-based FEM with ray directions d_c
    fprintf(['-'*ones(1,80) '\n']);
    fprintf('Step3: Ray-FEM, high frequency \n');
    
    tic;
    omega = high_omega(ti);
    wpml = high_wpml(ti);                % width of PML
    sigmaMax = 25/wpml;                 % Maximun absorbtion
    ray = numray1;
    
    % build linear system
    A = assemble_Helmholtz_matrix_RayFEM(mnode,melem,omega,wpml,sigmaMax,speed,ray,fquadorder);
    source = @(p) nodal_basis(xs,ys,h,p);
    b = assemble_RHS_RayFEM(mnode,melem,omega,wpml,sigmaMax,source,speed,ray,fquadorder);
    b = b/(h*h); %normalize b
    
    uh1 = RayFEM_direct_solver(mnode,melem,A,b,omega,ray,speed);
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
    cnumray = zeros(cN,Nray);
    
    
    % NMLA
    tic;
    for i = 1:cN
        x0 = cnode(i,1);  y0 = cnode(i,2);
        d0 = sqrt((x0-xs)^2 + (y0-ys)^2);
        c0 = speed(cnode(i,:));
        if d0 > r
            Rest = d0;
            cnumray_angle = NMLA(x0,y0,c0,omega,Rest,mnode,melem,uh1,ux,uy,[],1/5,Nray,'num',sec_opt,plt);
            cnumray(i) = exp(1i*cnumray_angle);
        else
            cnumray(i) = ex_ray([x0,y0],xs,ys,1);
        end
    end
    numray2 = interpolation(cnode,celem,node,cnumray);
    toc;
    
    % compute the ray errors
    exray = ex_ray(node,xs,ys,1);
    xx = node(:,1) - xs;   yy = node(:,2) - ys;
    rr = sqrt(xx.*xx + yy.*yy);
    numray2(rr<=r) = exray(rr<=r);
    rayerr2 = numray2 - exray;
    
    
    %% Step 5: Solve the original Helmholtz equation by Ray-based FEM with ray directions d_o
    fprintf(['-'*ones(1,80) '\n']);
    fprintf('Step5: Ray-FEM, high frequency \n');
    
    tic;
    omega = high_omega(ti);
    wpml = high_wpml(ti);                % width of PML
    sigmaMax = 25/wpml;                 % Maximun absorbtion
    ray = numray2;
    
    % build linear system
    A = assemble_Helmholtz_matrix_RayFEM(node,elem,omega,wpml,sigmaMax,speed,ray,fquadorder);
    source = @(p) nodal_basis(xs,ys,h,p);
    b = assemble_RHS_RayFEM(node,elem,omega,wpml,sigmaMax,source,speed,ray,fquadorder);
    b = b/(h*h); %normalize b
    
    uh = RayFEM_direct_solver(node,elem,A,b,omega,ray,speed);   
    toc;
    
    
    %% SFEM
    A = assemble_Helmholtz_matrix_SFEM(node,elem,omega,wpml,sigmaMax,speed,fquadorder);
    b = assemble_RHS_SFEM(node,elem, @(x)nodal_basis(xs,ys,h,x),fquadorder);
    b = b/(h*h/2);
    [~,~,isBdNode] = findboundary(elem);
    u_std = zeros(N,1);
    u_std(freeNode) = A(freeNode,freeNode)\b(freeNode);
    
end

totaltime = toc(tstart);
fprintf('\n\nTotal running time: % d minutes \n', totaltime/60);


%% map to polar
figure(30);
r1 = 0.811;
r2 = 0.811 + 2.1*high_wl;
theta1 = pi/4 - pi/8;
theta2 = pi/4 + pi/8;
subplot(2,1,1);
mapto_polar(node,elem,omega,speed,v,ray,xs,ys,r1,r2,theta1,theta2);

subplot(2,1,2);
mapto_polar(node,elem,omega,speed,u_std,0*ray,xs,ys,r1,r2,theta1,theta2);
