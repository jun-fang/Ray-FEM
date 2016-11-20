%% Test One Point Source Problem (outside domain): Iterative Idea
% uexact = sqrt(k)*besselh(0,1,k*sqrt((x-2)^2 + (y-2)^2))
clear;
fileID = fopen('result2.txt','a');


%% Load source data
global xs ys omega a;
xs = 20; ys = 20;
pde = Helmholtz_data_2;
fprintf(fileID,'\n\nFour point sources problem: \n\n');
fprintf(fileID,'(%d,%d),(%d,%d),(%d,%d),(%d,%d) \n\n', xs,ys,-xs,ys,xs,-ys,-xs,-ys);


%% Set up
plt = 0;                   % show solution or not
fquadorder = 9;            % numerical quadrature order
solver = 'DIR';            % linear system solver
pct = 1/8;
Nray = 4;
sec_opt = 0;               % NMLA second order correction or not

rec_N = 6;                 % we test rec_N examples

% record h and omega
rec_omega = zeros(1,rec_N);
rec_h = rec_omega;

% record the error of S-FEM for low frequency equation
rec_low_err = rec_omega;

% record the L2 error of the numerical angle estimation
rec_ang_err1 = rec_omega;
rec_ang_err2 = rec_omega;

% record the error of Numerical Ray-based FEM
rec_NR_err1 = rec_omega;
rec_NR_err2 = rec_omega;

% record the error of Exact Ray-based FEM
rec_ER_err = rec_omega;


% lg_a,md_a,sm_a have to be chosen to match the mesh size h and ch
% (for example md_a/ch = integer), and the real radius in NMLA as well
% (r1+md_a < lg_a, r2 + sm_a < md_a)
% lg_a = 2.1;  md_a = 1.1;  sm_a = 0.5;
% Rest1 = 20;  Rest2 = 20;
NPW = 6;

% cp_omega = [10 20 40 60 80 120 160]*pi;
% cp_omega = [10 20 40 80]*pi;
cp_omega = [10 20 30 40 60 80 120 160]*pi;

Rest1 = 30; Rest2 = 30;
sm_a = 1/2;
low_r = NMLA_radius(sqrt(cp_omega),Rest1);
% low_r(1) = NMLA_radius(sqrt(cp_omega(1)),Rest1);
high_r = NMLA_radius(cp_omega,Rest2);
M_a = sm_a + ceil(high_r/0.1)*0.1;
L_a = sm_a + M_a + ceil(low_r/0.1)*0.1;
fprintf('L_a =   ');
fprintf( '&  %1.2f  ',L_a);
fprintf('\nM_a =   ');
fprintf( '&  %1.2f  ',M_a);


tstart = tic;
for rec_i = 1: rec_N
    high_omega = cp_omega(rec_i);
    low_omega = sqrt(high_omega);
    h = 1/(NPW*round(high_omega/(2*pi)));
    ch = 1/(10*round(low_omega/(2*pi)));
    
    lg_a = L_a(rec_i);  md_a = M_a(rec_i);
%     if high_omega >= 20*pi
%         Rest1 = 10; Rest2 = 10;
%     end
    
%     if (high_omega >= 10*pi) && (high_omega < 20*pi)
%         lg_a = 1.9;  md_a = 1.1;  sm_a = 0.5;
%         Rest1 = 12;  Rest2 = 10;
%     elseif (high_omega >= 20*pi) && (high_omega < 40*pi)
%         lg_a = 1.5;  md_a = 0.8;  sm_a = 0.5;
%         Rest1 = 10;  Rest2 = 10;
%     elseif (high_omega >= 40*pi) && (high_omega < 60*pi)
%         lg_a = 1.3;  md_a = 0.8;  sm_a = 0.5;
%         Rest1 = 10;  Rest2 = 10;
%     elseif (high_omega >= 60*pi) && (high_omega < 80*pi)
%         lg_a = 1.2;  md_a = 0.7;  sm_a = 0.5;
%         Rest1 = 10;  Rest2 = 10;
%     elseif (high_omega >= 80*pi) && (high_omega < 120*pi)
%         lg_a = 1.2;  md_a = 0.7;  sm_a = 0.5;
%         Rest1 = 10;  Rest2 = 10;
%     elseif (high_omega >= 120*pi) 
%         lg_a = 1.1;  md_a = 0.7;  sm_a = 0.5;
%         Rest1 = 10;  Rest2 = 10;
%     end
    
%     if (high_omega >= 20*pi -eps) && (high_omega < 40*pi -eps)
%         lg_a = 1.9;  md_a = 1.0;  sm_a = 0.5;
%     elseif (high_omega >= 40*pi -eps) && (high_omega < 80*pi -eps)
%         lg_a = 1.5;  md_a = 0.8;  sm_a = 0.5;
%     elseif (high_omega >= 80*pi -eps) && (high_omega < 160*pi -eps)
%         lg_a = 1.4;  md_a = 0.8;  sm_a = 0.5;
%     elseif (high_omega >= 160*pi -eps) 
%         lg_a = 1.2;  md_a = 0.7;  sm_a = 0.5;
%     end
    
    fprintf(['-'*ones(1,80) '\n']);
    fprintf(['-'*ones(1,80) '\n']);
    fprintf('\ncase %d: \nomega/(2*pi) = %d,   1/h = %d   1/ch = %d,  NPW = %d \n\n',...
        rec_i, round(high_omega/(2*pi)), round(1/h), round(1/ch), NPW);
    fprintf('lg_a = %.2f,   md_a = %.2f,   sm_a = %.2f,   Rest1 = %.2f,  Rest2 = %.2f\n',...
        lg_a, md_a, sm_a, Rest1, Rest2);
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
    [u_std,~,~,rec_low_err(rec_i)] = Standard_FEM_IBC(lnode,lelem,omega,pde,fquadorder,solver,plt);
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
    Rest = Rest1;
    for i = 1:cN
        x0 = cnode(i,1);  y0 = cnode(i,2);
        c0 = pde.speed(cnode(i,:));
        [cnumray(i,:)] = NMLA_2D_correction(x0,y0,c0,omega,Rest,lnode,lelem,u_std,ux,uy,pde,pct,Nray,'num',sec_opt,plt);
    end
    numray1 = interpolation(cnode,celem,mnode,cnumray);
    toc;
    
    ray = pde.ray(mnode);
    ray = ray(:);
    ray = [real(ray), imag(ray)];
    ray_dir = atan2(ray(:,2),ray(:,1));
    neg_index = find(ray_dir<0);
    ray_dir(neg_index) = ray_dir(neg_index) + 2*pi;
    
    diffang = numray1(:) - ray_dir;
    rec_ang_err1(rec_i) = h*norm(diffang,2);
    
    numray = exp(1i*numray1);
    
    
    if(1)
    %% Step 3: Solve the original Helmholtz equation by Ray-based FEM with ray directions d_c
    fprintf(['\n' '-'*ones(1,80) '\n']);
    fprintf('\nStep3 \n');
    tic;
    omega = high_omega;
    uh1 = Ray_FEM_IBC_1(mnode,melem,omega,pde,numray,fquadorder,plt);
    exu1 = pde.ex_u(mnode);
    du1 = uh1-exu1;
    rec_NR_err1(rec_i) = h*norm(du1,2);
    toc;
%     uh1 = pde.ex_u(mnode);
    
    idx = find((abs(mnode(:,1))<=sm_a+eps).*(abs(mnode(:,2))<=sm_a+eps));
    rec_NR_err1(rec_i) = h*norm(du1(idx),2);
    dang = reshape(diffang,mN,Nray);
    dang = dang(idx,:);
    rec_ang_err1(rec_i) = h*norm(dang(:),2);

    
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
    Rest = Rest2;
    for i = 1:cN
        x0 = cnode(i,1);  y0 = cnode(i,2);
        c0 = pde.speed(cnode(i,:));
        [cnumray(i,:)] = NMLA_2D_correction(x0,y0,c0,omega,Rest,mnode,melem,uh1,ux,uy,pde,pct,Nray,'num',sec_opt,plt);
    end
    numray2 = interpolation(cnode,celem,node,cnumray);
    toc;
    
    ray = pde.ray(node);
    ray = ray(:);
    ray = [real(ray), imag(ray)];
    ray_dir = atan2(ray(:,2),ray(:,1));
    neg_index = find(ray_dir<0);
    ray_dir(neg_index) = ray_dir(neg_index) + 2*pi;
    
    diffang = numray2(:) - ray_dir;
    rec_ang_err2(rec_i) = h*norm(diffang,2);
    
    numray = exp(1i*numray2);
    
    end
    
    
    if (1)
    %% Step 5: Solve the original Helmholtz equation by Ray-based FEM with ray directions d_o
    fprintf(['\n' '-'*ones(1,80) '\n']);
    fprintf('\nStep5 \n');
    tic;
    omega = high_omega;
    uh2 = Ray_FEM_IBC_1(node,elem,omega,pde,numray,fquadorder,plt);
    exu2 = pde.ex_u(node);
    rec_NR_err2(rec_i) = h*norm(uh2-exu2,2);
    toc;
    
    end
    
    
    %% Exact Ray-based FEM:
    if (0)
        fprintf('\nExact Ray-based FEM: \n');
        tic;
        ray = pde.ray(node);
        uh = Ray_FEM_IBC_1(node,elem,omega,pde,ray,fquadorder,plt);
        exu = pde.ex_u(node);
        rec_ER_err(rec_i) = h*norm(uh-exu,2);
        toc;
    end
    
    
end

totaltime = toc(tstart);
fprintf('\n\nTotal running time: % d minutes \n', totaltime/60);



%% record and print results
rec = [rec_omega;rec_h;rec_low_err;rec_ER_err;rec_ang_err1;rec_ang_err2;rec_NR_err1;rec_NR_err2];

save('result2.mat','rec_omega','rec_h','rec_low_err','rec_ang_err1','rec_ang_err2',...
    'rec_NR_err1','rec_NR_err2','rec_ER_err');


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
fprintf( fileID,'\nS-FEM for low frequency equation:\n\n');
fprintf( fileID,'\n\nRelative L2 error:      ');
fprintf( fileID,'&  %1.2d  ',rec_low_err);

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


fprintf( fileID,['\n' '-'*ones(1,80) '\n']);
fprintf( fileID,'\nExact Ray-based FEM:\n\n');
fprintf( fileID,'\n\nRelative L2 error:      ');
fprintf( fileID,'&  %1.2d  ',rec_ER_err);



fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('S-FEM for low frequency equation:\n');
fprintf('\nRelative L2 error:      ');
fprintf( '&  %1.2d  ',rec_low_err);

fprintf( ['\n\n' '-'*ones(1,80) '\n']);
fprintf( 'Numerical Ray-based FEM:\n\n');
fprintf( 'Angle L2 error 1:       ');
fprintf( '&  %1.2d  ',rec_ang_err1);
fprintf( '\n\nAngle L2 error 2:       ');
fprintf( '&  %1.2d  ',rec_ang_err2);
fprintf( '\n\nRelative L2 error 1:    ');
fprintf( '&  %1.2d  ',rec_NR_err1);
fprintf( '\n\nRelative L2 error 2:    ');
fprintf( '&  %1.2d  ',rec_NR_err2);

fprintf( ['\n\n' '-'*ones(1,80) '\n']);
fprintf( 'Exact Ray-based FEM:');
fprintf( '\n\nRelative L2 error:      ');
fprintf( '&  %1.2d  ',rec_ER_err);
fprintf( ['\n' '-'*ones(1,80) '\n']);


