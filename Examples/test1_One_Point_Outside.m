%% Convergence test for homogenenous case with exact ray information

clear;
addpath(genpath('../../ifem/'));
addpath('../Methods/');
addpath('../Cutoff_Functions/')
addpath('../Plots_Prints/');


%% Set up
xs = 0; ys = 0;                     % source location 
epsilon = 50/(80*pi);               % cut-off parameter   
speed = @(p) ones(size(p(:,1)));    % wave speed

wpml = 0.1;                         % width of PML  
sigmaMax = 25/wpml;                 % absorption of PML  
fquadorder = 3;                     % numerical quadrature order 
a = 1/2;                            % computational domain [-a,a]^2

nt = 4;                             % number of tests
errors = zeros(1,nt);                  % record errors
omegas = pi*[120,160,240,320];      % omega's
NPW = 4;                            % grid number per wavelength


%% Tests
for ii = 1:nt
    
    omega = omegas(ii);
    wl = 2*pi/omega;
    h = 1/round(1/(wl/NPW));
    [node,elem] = squaremesh([-a,a,-a,a],h);
    fprintf('Case %d: omega/pi = %d, NPW = %d, 1/h = %d\n', ii, round(omega/pi), NPW, round(1/h));
    
    %% Exact ray information
    x = node(:,1);  y = node(:,2);
    rr = sqrt((x-xs).^2 + (y-ys).^2);
    ray = atan2(y-ys, x-xs);
    ray = exp(1i*ray).*(rr>10*eps);
    
    %% Assemble and solve the system Au = b
    option = 'homogeneous';
    A = assemble_Helmholtz_matrix_RayFEM(node,elem,omega,wpml,sigmaMax,speed,ray,fquadorder);
    b = assemble_RHS_RayFEM_with_ST(node,elem,xs,ys,omega,epsilon,wpml,sigmaMax,ray,speed,fquadorder,option);
    u = direct_solver(node,elem,A,b,omega,ray,speed);
       
    %% Get the exact solution
    ub = 1i/4*besselh(0,1,omega*rr);
    cf = cutoff(epsilon,2*epsilon,node,xs,ys);
    uex = (1-cf).*ub;
    uex(rr<epsilon) = 0;
    
    %% Compute relative L2 error 
    du = u - uex;
    idx = find( ~( (x<=max(x)-wpml).*(x>= min(x)+wpml)...
        .*(y<= max(y)-wpml).*(y>= min(y)+wpml) ) ); % index on PML
    du(idx) = 0;  uex(idx) = 0;
    errors(ii) = norm(du)/norm(uex);
    toc;
end


%% plot
figure(22);
show_convergence_rate(omegas(1:nt),errors,'omega','||u - u_h||_{L^2(\Omega)}/||u||_{L^2(\Omega)}');

























































%% Test One Point Source Problem (outside domain): Iterative Idea
% uexact = sqrt(k)*besselh(0,1,k*sqrt((x-2)^2 + (y-2)^2))
clear;
fileID = fopen('result1_iter.txt','a');


%% Load source data
pde = Helmholtz_data1;
fprintf(fileID,'\n\n One point source problem (outside domain): \n\n');
fprintf(fileID,'u_ex = sqrt(k)*besselh(0,1,k*sqrt((x-2)^2 + (y-2)^2)) \n\n');
xs = 2; ys = 2;            % point source location


%% Set up
plt = 0;                   % show solution or not
fquadorder = 9;            % numerical quadrature order
solver = 'DIR';            % linear system solver
Nray = 1;                  % one ray direction
sec_opt = 0;               % NMLA second order correction or not

rec_N = 6;                 % we test rec_N examples

% record h and omega
rec_omega = zeros(1,rec_N);
rec_h = rec_omega;

% record the error and condition number of Standard FEM
rec_S_err = rec_omega;
rec_S_cond = rec_omega;

% record the L2 error of the numerical angle estimation
rec_ang_err1 = rec_omega;
rec_ang_err2 = rec_omega;

% record the error and condition number of Numerical Ray-based FEM
rec_NR_err1 = rec_omega;
rec_NR_err2 = rec_omega;
rec_NR_cond1 = rec_omega;
rec_NR_cond2 = rec_omega;

% record the error and condition number of Exact Ray-based FEM
rec_ER_err = rec_omega;
rec_ER_cond = rec_omega;

% record the interpolation error with exact ray information
rec_int_err = rec_omega;

global omega;
global a;
% lg_a,md_a,sm_a have to be chosen to match the mesh size h and ch
% (for example md_a/ch = integer), and the real radius in NMLA as well
% (r1+md_a < lg_a, r2 + sm_a < md_a)
lg_a = 5/4;
md_a = 3/4;
sm_a = 1/2;
Rest = 1;
NPW = 6;

% cp_omega = [20 40 80 160]*2*pi;
cp_omega = [20 40 60 80 120 160 200 260 320 400]*pi;
tstart = tic;
for rec_i = 1: rec_N
    high_omega = cp_omega(rec_i);
    low_omega = sqrt(high_omega);
    h = 1/(NPW*round(high_omega/(2*pi)));
    ch = 1/(2*NPW*round(low_omega/(2*pi)));
    
    if (high_omega >= 100*pi) && (high_omega < 160*pi)
        lg_a = 1;
        md_a = 2/3;
        sm_a = 1/2;
    elseif high_omega >= 160*pi
        lg_a = 7/8;
        md_a = 5/8;
        sm_a = 1/2;
    end
    
    fprintf(['-'*ones(1,80) '\n']);
    fprintf(['-'*ones(1,80) '\n']);
    fprintf('\ncase %d: \nomega/(2*pi) = %d,   1/h = %d   1/ch = %d,  NPW = %d \n\n',...
        rec_i, round(high_omega/(2*pi)), 1/h, 1/ch, NPW);
    rec_omega(rec_i) = high_omega;
    rec_h(rec_i) = h;
    
    
    %% Step 1: Solve the Hemholtz equation with the same source but with a relative low frequency sqrt(\omega) by Standard FEM, mesh size \omega*h = constant
    fprintf(['\n' '-'*ones(1,80) '\n']);
    fprintf('Numerical Ray-based FEM: \n\n');
    fprintf('Step1 \n');
    tic;
    omega = low_omega;
    a = lg_a;
    [lnode,lelem] = squaremesh([-a,a,-a,a],h);
    [u_std,~,~,~] = Standard_FEM_IBC(lnode,lelem,omega,pde,fquadorder,solver,plt);
    toc;
    
    
    %% Step 2: Use NMLA to find ray directions d_c with low frequency sqrt(\omega)
    fprintf(['\n' '-'*ones(1,80) '\n']);
    fprintf('\nStep2 \n');
    
    lN = size(lnode,1);
    ln = round(sqrt(lN));
    n = ln;
    
    uu = reshape(u_std,n,n);
    uux = uu;   uuy = uu;   % numerical Du
    
    uux(:,2:n-1) = (uu(:,3:n) - uu(:,1:n-2))/(2*h);
    uux(:,n) = 2*uux(:,n-1) - uux(:,n-2);
    uux(:,1) = 2*uux(:,2) - uux(:,3);
    ux = uux(:);
    
    uuy(2:n-1,:) = (uu(3:n,:) - uu(1:n-2,:))/(2*h);
    uuy(n,:) = 2*uuy(n-1,:) - uuy(n-2,:);
    uuy(1,:) = 2*uuy(2,:) - uuy(3,:);
    uy = uuy(:);
    
    a = md_a;
    [mnode,melem] = squaremesh([-a,a,-a,a],h);
    mN = size(mnode,1);
    mn = round(sqrt(mN));
    
    [cnode,celem] = squaremesh([-a,a,-a,a],ch);
    cN = size(cnode,1);
    cnumray = zeros(cN,Nray);
    cr = zeros(cN,Nray);
    
    fprintf('NMLA time: \n');
    tic;
    for i = 1:cN
        x0 = cnode(i,1);  y0 = cnode(i,2);
        c0 = pde.speed(cnode(i,:));
        [cnumray(i)] = NMLA_2D_2nd(x0,y0,c0,omega,Rest,lnode,lelem,u_std,ux,uy,pde,1/5,Nray,'num',sec_opt,plt);
    end
    numray1 = interpolation(cnode,celem,mnode,cnumray);
    toc;
    
    exray = ex_ray_angle(mnode,xs,ys);
    diffang1 = numray1 - exray;
    rec_ang_err1(rec_i) = h*norm(diffang1,2)/(h*norm(exray,2));
    numray = exp(1i*numray1);
    
    
    %% Step 3: Solve the original Helmholtz equation by Ray-based FEM with ray directions d_c
    fprintf(['\n' '-'*ones(1,80) '\n']);
    fprintf('\nStep3 \n');
    tic;
    omega = high_omega;
    [uh1,A,~,~,rel_L2_err] = Ray_FEM_IBC_1(mnode,melem,omega,pde,numray,fquadorder,plt);
    toc;
    rec_NR_err1(rec_i) = rel_L2_err;
    %     rec_NR_cond1(rec_i) = condest(A);
    
    
    %% Step 4: NMLA to find original ray directions d_o with wavenumber k
    fprintf(['\n' '-'*ones(1,80) '\n']);
    fprintf('\nStep4 \n');
    
    a = sm_a;
    [node,elem] = squaremesh([-a,a,-a,a],h);
    
    [cnode,celem] = squaremesh([-a,a,-a,a],ch);
    cN = size(cnode,1);
    cnumray = zeros(cN,Nray);
    
    n = mn;
    uu = reshape(uh1,n,n);
    uux = uu;   uuy = uu;   % numerical Du
    
    uux(:,2:n-1) = (uu(:,3:n) - uu(:,1:n-2))/(2*h);
    uux(:,n) = 2*uux(:,n-1) - uux(:,n-2);
    uux(:,1) = 2*uux(:,2) - uux(:,3);
    ux = uux(:);
    
    uuy(2:n-1,:) = (uu(3:n,:) - uu(1:n-2,:))/(2*h);
    uuy(n,:) = 2*uuy(n-1,:) - uuy(n-2,:);
    uuy(1,:) = 2*uuy(2,:) - uuy(3,:);
    uy = uuy(:);
    
    fprintf('NMLA time: \n');
    tic;
    for i = 1:cN
        x0 = cnode(i,1);  y0 = cnode(i,2);
        c0 = pde.speed(cnode(i,:));
        [cnumray(i,:)] = NMLA_2D_2nd(x0,y0,c0,omega,Rest,mnode,melem,uh1,ux,uy,pde,1/5,Nray,'num',sec_opt,plt);
    end
    numray2 = interpolation(cnode,celem,node,cnumray);
    toc;
    
    exray = ex_ray_angle(node,xs,ys);
    diffang2 = numray2 - exray;
    rec_ang_err2(rec_i) = h*norm(diffang2,2)/(h*norm(exray,2));
    numray = exp(1i*numray2);
    
    
    %% Step 5: Solve the original Helmholtz equation by Ray-based FEM with ray directions d_o
    fprintf(['\n' '-'*ones(1,80) '\n']);
    fprintf('\nStep5 \n');
    tic;
    omega = high_omega;
    [uh2,A,~,~,rel_L2_err] = Ray_FEM_IBC_1(node,elem,omega,pde,numray,fquadorder,plt);
    toc;
    rec_NR_err2(rec_i) = rel_L2_err;
    %     rec_NR_cond2(rec_i) = condest(A);
    
    
    %% Standard FEM
    if (0)
        fprintf('\nStandard FEM: \n');
        [~,A,~,rel_L2_err] = Standard_FEM_IBC(node,elem,omega,pde,fquadorder,solver,plt);
        rec_S_err(rec_i) = rel_L2_err;
        rec_S_cond(rec_i) = condest(A);
    end
    
    
    %% Exact Ray-based FEM:
    if (0)
        fprintf('\nExact Ray-based FEM: \n');
        tic;
        ray = pde.ray(node);
        [~,A,~,~,rel_L2_err] = Ray_FEM_IBC_1(node,elem,omega,pde,ray,fquadorder,plt);
        toc;
        rec_ER_err(rec_i) = rel_L2_err;
        %         rec_ER_cond(rec_i) = condest(A);
        
        coeff = pde.int_coe(node);
        c = pde.speed(node);
        [rec_int_err(rec_i)] = Ray_FEM_L2_Error(coeff,node,elem,omega,c,pde.ex_u,ray,fquadorder);
    end
    
    
end

totaltime = toc(tstart);
fprintf('\n\nTotal running time: % d minutes \n', totaltime/60);






%% record and print results
rec = [rec_omega;rec_h;rec_S_err;rec_S_cond;...
    rec_ang_err1;rec_ang_err2;rec_NR_err1;rec_NR_err2;rec_NR_cond1;rec_NR_cond2;...
    rec_ER_cond;rec_ER_err;rec_int_err];

save('result1.mat','rec_omega','rec_h','rec_S_err','rec_S_cond','rec_ang_err1',...
    'rec_ang_err2','rec_NR_err1','rec_NR_err2','rec_NR_cond1','rec_NR_cond2',...
    'rec_ER_cond','rec_ER_err','rec_int_err');


fprintf( fileID,['\n' '-'*ones(1,80) '\n']);
fprintf( fileID,'omega:                  ');
fprintf( fileID,'&  %.2e  ',rec_omega );
fprintf( fileID,'\nomega/2pi:              ');
fprintf( fileID,'&  %.2e  ',rec_omega/(2*pi) );
fprintf( fileID,'\n\nGrid size h:            ');
fprintf( fileID,'&  %.2e  ',rec_h);
fprintf( fileID,'\n1/h:                    ');
fprintf( fileID,'&  %.2e  ',1./rec_h);

fprintf( fileID,['\n' '-'*ones(1,80) '\n']);
fprintf( fileID,'\nNumerical Ray-based FEM:\n\n');
fprintf( fileID,'Angle L2 error 1:       ');
fprintf( fileID,'&  %1.2d  ',rec_ang_err1);
fprintf( fileID,'\n\nAngle L2 error 2:       ');
fprintf( fileID,'&  %1.2d  ',rec_ang_err2);
fprintf( fileID,'\n\nRelative L2 error 1:    ');
fprintf( fileID,'&  %1.2d  ',rec_NR_err1);
fprintf( fileID,'\n\nRelative L2 error 2:    ');
fprintf( fileID,'&  %1.2d  ',rec_NR_err2);
fprintf( fileID,'\n\nCondition number 1:     ');
fprintf( fileID,'&  %1.2d  ',rec_NR_cond1);
fprintf( fileID,'\n\nCondition number 2:     ');
fprintf( fileID,'&  %1.2d  ',rec_NR_cond2);

fprintf( fileID,['\n' '-'*ones(1,80) '\n']);
fprintf( fileID,'\nExact Ray-based FEM:\n\n');
fprintf( fileID,'Condition number:       ');
fprintf( fileID,'&  %1.2d  ',rec_ER_cond);
fprintf( fileID,'\n\nRelative L2 error:      ');
fprintf( fileID,'&  %1.2d  ',rec_ER_err);

fprintf( fileID,['\n\n' '-'*ones(1,80) '\n']);
fprintf( fileID,'Interpolation error:    ');
fprintf( fileID,'&  %1.2d  ',rec_int_err);

fprintf( fileID,['\n' '-'*ones(1,80) '\n']);
fprintf( fileID,'\nStandard FEM:\n\n');
fprintf( fileID,'Condition number:       ');
fprintf( fileID,'&  %1.2d  ',rec_S_cond);
fprintf( fileID,'\n\nRelative L2 error:      ');
fprintf( fileID,'&  %1.2d  ',rec_S_err);



fprintf( ['\n' '-'*ones(1,80) '\n']);
fprintf( '\nNumerical Ray-based FEM:\n\n');
fprintf( 'Angle L2 error 1:       ');
fprintf( '&  %1.2d  ',rec_ang_err1);
fprintf( '\n\nAngle L2 error 2:       ');
fprintf( '&  %1.2d  ',rec_ang_err2);
fprintf( '\n\nRelative L2 error 1:    ');
fprintf( '&  %1.2d  ',rec_NR_err1);
fprintf( '\n\nRelative L2 error 2:    ');
fprintf( '&  %1.2d  ',rec_NR_err2);
fprintf( ['\n' '-'*ones(1,80) '\n']);

fprintf( ['\n' '-'*ones(1,80) '\n']);
fprintf( '\nExact Ray-based FEM:');
fprintf( '\n\nRelative L2 error:      ');
fprintf( '&  %1.2d  ',rec_ER_err);
fprintf( ['\n' '-'*ones(1,80) '\n']);