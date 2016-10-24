%% Point source problem in homogeneous medium
% uexact = sqrt(k)*besselh(0,1,k*sqrt((x-xs)^2 + (y-ys)^2))
clear;
fileID = fopen('test1_result.txt','a');


%% Load source data
pde = Helmholtz_data_point_source;
global omega xs ys a;
xs = 0; ys = 0;                     % point source location
speed = @(x) ones(size(x,1),1);     % medium speed

fprintf(fileID,['-'*ones(1,80) '\n']);
fprintf(fileID,'Point source problem in homogeneous medium:\n\n  ');
fprintf(fileID,'u_ex = sqrt(k)*besselh(0,1,k*sqrt((x-%d)^2 + (y-%d)^2)) \n\n',xs,ys);




%% Set up
plt = 0;                   % show solution or not
fquadorder = 3;            % numerical quadrature order
solver = 'DIR';            % linear system solver
Nray = 1;                  % one ray direction
sec_opt = 1;               % NMLA second order correction or not
epsilon = 0.137;

NPW = 6;                   % number of points per wavelength
test_num = 6;              % we test test_num examples

% frequency
high_omega = [40 80 120 160 240 320]*pi;
low_omega = sqrt(high_omega);

% error
low_rayerr = 0*high_omega;
high_rayerr = 0*high_omega;
l2_err = 0*high_omega;
rel_l2_err = 0*high_omega;
rel_max_err = 0*high_omega;
max_err = 0*high_omega;

% wavelength
high_wl = 2*pi./high_omega;
low_wl = 2*pi./low_omega;

% mesh size
fh = 1./(NPW*round(high_omega/(2*pi)));      % fine mesh size
ch = 1./(20*round(low_omega/(2*pi)));        % coarse mesh size

% width of PML
% high_wpml = ch*ceil(high_wl/ch);
high_wpml = 0.08*ones(size(high_omega)); %fh.*ceil(high_wl./fh);
low_wpml = fh.*ceil(low_wl./fh);


%% Generate the domain sizes
sd = 1/2;
Rest = sqrt(2)*sd;

high_r = NMLA_radius(high_omega,Rest);
md = sd + high_r + high_wpml;
md = ceil(md*10)/10;

Rest = sqrt(2)*md;
low_r = NMLA_radius(low_omega,Rest);
ld = md + low_r + low_wpml;
ld = ceil(ld*10)/10;


%% Test
tstart = tic;
for ti = 1: test_num
    omega = high_omega(ti);
    h = fh(ti);  h_c = ch(ti);
    fprintf(['-'*ones(1,80) '\n']);
    fprintf('\ncase %d: \nomega/(2*pi) = %d,   1/h = %d   1/h_c = %d,  NPW = %d \n',...
        ti, round(omega/(2*pi)), 1/h, 1/h_c, NPW);
    
    
    %% Step 1: Solve the Hemholtz equation with the same source but with a relative low frequency sqrt(\omega) by Standard FEM, mesh size \omega*h = constant
    fprintf(['\n' '-'*ones(1,80) '\n']);
    fprintf('Step1: S-FEM, low frequency \n');
    tic;
    omega = low_omega(ti);
    a = ld(ti);
    wpml = low_wpml(ti);                % width of PML
    sigmaMax = 25/wpml;                 % Maximun absorbtion
    [lnode,lelem] = squaremesh([-a,a,-a,a],h);
    [u_std] = Standard_FEM_PML_PointSource(lnode,lelem,omega,wpml,sigmaMax,xs,ys,speed,fquadorder,plt);
    toc;
    
    
    %% Step 2: Use NMLA to find ray directions d_c with low frequency sqrt(\omega)
    fprintf(['\n' '-'*ones(1,80) '\n']);
    fprintf('\nStep2: NMLA, low frequency \n');
    
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
    
    a = md(ti);
    [mnode,melem] = squaremesh([-a,a,-a,a],h);
    mN = size(mnode,1);  mn = round(sqrt(mN));
    
    [cnode,celem] = squaremesh([-a,a,-a,a],h_c);
    cN = size(cnode,1);
    cnumray_angle = zeros(cN,Nray);
    
    tic;
    for i = 1:cN
        x0 = cnode(i,1);  y0 = cnode(i,2);
        r0 = sqrt((x0-xs)^2 + (y0-ys)^2);
        c0 = pde.speed(cnode(i,:));
        if r0>2*epsilon
            Rest = r0;
            [cnumray_angle(i)] = NMLA_2D_2nd(x0,y0,c0,omega,Rest,lnode,lelem,u_std,ux,uy,pde,1/5,Nray,'num',sec_opt,plt);
        else
            cnumray_angle(i) = ex_ray_angle([x0,y0],xs,ys);
        end
    end
    cnumray = exp(1i*cnumray_angle);
    numray1 = interpolation(cnode,celem,mnode,cnumray);
    toc;
    
    exray = ex_ray(mnode,xs,ys,1);
    mr = sqrt((mnode(:,1)-xs).^2 + (mnode(:,2)-ys).^2);
    numray1 = numray1.*(mr>2*epsilon) + exray.*(mr<=2*epsilon);
    rayerr1 = numray1 - exray;
    low_rayerr(ti) = norm(rayerr1)*h/(norm(exray)*h);
    numray = numray1;
    
    
    %% Step 3: Solve the original Helmholtz equation by Ray-based FEM with ray directions d_c
    fprintf(['\n' '-'*ones(1,80) '\n']);
    fprintf('\nStep3: Ray-FEM, high frequency \n');
    %     tic;
    omega = high_omega(ti);
    wpml = high_wpml(ti);                % width of PML
    sigmaMax = 25/wpml;                 % Maximun absorbtion
    [uh1] = Ray_FEM_PML_1_PointSource(mnode,melem,omega,wpml,sigmaMax,xs,ys,speed,numray,fquadorder,plt);
    %     toc;
    
    
    
    %% Step 4: NMLA to find original ray directions d_o with wavenumber k
    fprintf(['\n' '-'*ones(1,80) '\n']);
    fprintf('\nStep4: NMLA, high frequency \n');
    
    a = sd;
    [node,elem] = squaremesh([-a,a,-a,a],h);
    
    [cnode,celem] = squaremesh([-a,a,-a,a],h_c);
    cN = size(cnode,1);
    cnumray_angle = zeros(cN,Nray);
    
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
    
    tic;
    for i = 1:cN
        x0 = cnode(i,1);  y0 = cnode(i,2);
        r0 = sqrt((x0-xs)^2 + (y0-ys)^2);
        c0 = pde.speed(cnode(i,:));
        if r0>2*epsilon
            Rest = r0;
            [cnumray_angle(i)] = NMLA_2D_2nd(x0,y0,c0,omega,Rest,mnode,melem,uh1,ux,uy,pde,1/5,Nray,'num',sec_opt,plt);
        else
            cnumray_angle(i) = ex_ray_angle([x0,y0],xs,ys);
        end
    end
    cnumray = exp(1i*cnumray_angle);
    numray2 = interpolation(cnode,celem,node,cnumray);
    toc;
    
    exray = ex_ray(node,xs,ys,1);
    sr = sqrt((node(:,1)-xs).^2 + (node(:,2)-ys).^2);
    numray2 = numray2.*(sr>2*epsilon) + exray.*(sr<=2*epsilon);
    rayerr2 = numray2 - exray;
    high_rayerr(ti) = norm(rayerr2)*h/(norm(exray)*h);
    numray = numray2;
    
    
    
    %% Step 5: Solve the original Helmholtz equation by Ray-based FEM with ray directions d_o
    fprintf(['\n' '-'*ones(1,80) '\n']);
    fprintf('\nStep5: Ray-FEM, high frequency \n');
    tic;
    
    % Assembling
    wpml = 0.08;                % width of PML
    sigmaMax = 25/wpml;                 % Maximun absorbtion
    ray = numray;
    [A] = assemble_Helmholtz_matrix_with_ray_1(node,elem,omega,wpml,sigmaMax,speed,ray,fquadorder);
    b = assemble_RHS_with_ray_sing(node,elem,epsilon,omega,speed,ray,fquadorder);
    
    % Boundaries
    [~,~,isBdNode] = findboundary(elem);
    rep_isBdNode = repmat(isBdNode,1,Nray);
    isBdNode = rep_isBdNode(:);
    freeNode = find(~isBdNode);
    
    % Solve the linear system Au = b
    N = size(node,1); Ndof = N;
    v = zeros(Ndof,1);
    v(freeNode) = A(freeNode,freeNode)\b(freeNode);
    
    % Compute solution values at grid nodes
    grad = ray(:);
    grad = [real(grad),imag(grad)];
    repnode = repmat(node,Nray,1);
    temp = grad(:,1).*repnode(:,1) + grad(:,2).*repnode(:,2);
    
    k = omega./speed(node);           % wavenumber
    kk = repmat(k,1,Nray);
    u = v.*exp(1i*kk(:).*temp);
    u = reshape(u,N,Nray);
    u = sum(u,2);
    
    
    % Construct the solution
    xx = node(:,1)-xs;  yy = node(:,2)-ys;
    rr = sqrt(xx.^2 + yy.^2);
    ub = 1i/4*besselh(0,1,omega*rr);
    x_eps = cutoff(epsilon,2*epsilon,node);
    v = ub.*x_eps;
    us = u;
    uu = us + v;
    toc;
    
    % compute the L2 error
    du = ub -uu;
    x = node(:,1); y = node(:,2);
    rxy = sqrt((x-xs).^2 + (y-ys).^2);
    du(rxy<h/2) = 0;
    du(x>=a-wpml)=0; du(x<= -a+wpml) = 0;
    du(y>=a-wpml)=0; du(y<= -a+wpml) = 0;
    
    
    %% Compute the errors of wavefield at y=0.4
    px = 0.4;
    px = round(px/h)*h;
    indx = 1 + round((px+a)/h);
    n = round(sqrt(N));
    du = reshape(du,n,n);
    uex = reshape(ub,n,n);
    
    du = du(indx,:);
    eu = uex(indx,:);
    
    max_err(ti) = norm(du,inf);
    rel_max_err(ti) = norm(du,inf)/norm(eu,inf);
    l2_err(ti) = norm(du)*h;
    rel_l2_err(ti) = norm(du)/norm(eu);
    
    
end

totaltime = toc(tstart);
fprintf('\n\nTotal running time: % d minutes \n', totaltime/60);

figure(1);
subplot(1,2,1);
FJ_showrate(high_omega(1:test_num),low_rayerr(1:test_num));
subplot(1,2,2);
FJ_showrate(high_omega(1:test_num),high_rayerr(1:test_num));

figure(2);
subplot(2,2,1);
FJ_showrate(high_omega(1:test_num),max_err(1:test_num));
subplot(2,2,2);
FJ_showrate(high_omega(1:test_num),rel_max_err(1:test_num));
subplot(2,2,3);
FJ_showrate(high_omega(1:test_num),l2_err(1:test_num));
subplot(2,2,4);
FJ_showrate(high_omega(1:test_num),rel_l2_err(1:test_num));




%% print result
fprintf( fileID,['\n' '-'*ones(1,80) '\n']);
fprintf( fileID,'omega:                  ');
fprintf( fileID,'&  %.2e  ',high_omega );
fprintf( fileID,'\nomega/2pi:              ');
fprintf( fileID,'&  %.2e  ',high_omega/(2*pi) );
fprintf( fileID,'\n\nGrid size h:            ');
fprintf( fileID,'&  %.2e  ',fh);
fprintf( fileID,'\n1/h:                    ');
fprintf( fileID,'&  %.2e  ',1./fh);

fprintf( fileID,['\n' '-'*ones(1,80) '\n']);
fprintf( fileID,'Low freq ray error:     ');
fprintf( fileID,'&  %1.2d  ',low_rayerr);
fprintf( fileID,'\n\nHigh freq ray error:    ');
fprintf( fileID,'&  %1.2d  ',high_rayerr);
fprintf( fileID,'\n\nMax error:              ');
fprintf( fileID,'&  %1.2d  ',max_err);
fprintf( fileID,'\n\nRelative max error:     ');
fprintf( fileID,'&  %1.2d  ',rel_max_err);
fprintf( fileID,'\n\nL2 error:               ');
fprintf( fileID,'&  %1.2d  ',l2_err);
fprintf( fileID,'\n\nRelative L2 error:      ');
fprintf( fileID,'&  %1.2d  ',rel_l2_err);


fprintf( ['\n' '-'*ones(1,80) '\n']);
fprintf( 'omega:                  ');
fprintf( '&  %.2e  ',high_omega );
fprintf( '\nomega/2pi:              ');
fprintf( '&  %.2e  ',high_omega/(2*pi) );
fprintf( '\n\nGrid size h:            ');
fprintf( '&  %.2e  ',fh);
fprintf( '\n1/h:                    ');
fprintf( '&  %.2e  ',1./fh);

fprintf( ['\n' '-'*ones(1,80) '\n']);
fprintf( 'Low freq ray error:     ');
fprintf( '&  %1.2d  ',low_rayerr);
fprintf( '\n\nHigh freq ray error:    ');
fprintf( '&  %1.2d  ',high_rayerr);
fprintf( '\n\nMax error:              ');
fprintf( '&  %1.2d  ',max_err);
fprintf( '\n\nRelative max error:     ');
fprintf( '&  %1.2d  ',rel_max_err);
fprintf( '\n\nL2 error:               ');
fprintf( '&  %1.2d  ',l2_err);
fprintf( '\n\nRelative L2 error:      ');
fprintf( '&  %1.2d  ',rel_l2_err);


fprintf( ['\n' '-'*ones(1,80) '\n']);


