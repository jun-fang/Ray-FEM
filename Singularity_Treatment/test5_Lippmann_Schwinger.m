%% Test for Lippmann-Schwinger equation

clear;

addpath(genpath('../../ifem/'));
addpath('../Methods/');
addpath('../NMLA/');
addpath('../Cutoff_Functions/')
addpath('../Plots_Prints/');


%% Load source/wavespeed data
xs = -0.4; ys = -0.4;                     % point source location

sigma = 0.15;
xHet = 0.2;   yHet = 0.2;

nu = @(x,y) 0.9*exp( -1/(2*sigma^2)*((x-xHet).^2 + (y-yHet).^2) )...
    .*Lippmann_Schwinger_window(sqrt((x-xHet).^2 + (y-yHet).^2), 0.28,0.2  );

speed = @(p) 1./sqrt(1 + nu( p(:,1), p(:,2) ));    % wave speed
% speed = @(p) ones(size(p(:,1)));

%% Set up
plt = 0;                   % show solution or not
fquadorder = 3;            % numerical quadrature order
Nray = 1;                  % one ray direction
sec_opt = 0;               % NMLA second order correction or not
epsilon = 50/(80*pi);               % cut-off parameter


NPW = 4;                   % number of points per wavelength
test_num = 6;              % we test test_num examples

% frequency
high_omega = [120 160 240 320 480 640]*pi;
low_omega = 2*sqrt(high_omega);

% error
max_err = 0*high_omega;            % L_inf error of the numerical solution
rel_max_err = 0*high_omega;        % relative L_inf error of the numerical solution
l2_err = 0*high_omega;             % L_2 error of the numerical solution
rel_l2_err = 0*high_omega;         % relative L_2 error of the numerical solution


% wavelength
high_wl = 2*pi./high_omega;
low_wl = 2*pi./low_omega;

% mesh size
fh = 1./(NPW*round(high_omega/(2*pi)));      % fine mesh size
ch = 1./(20*round(low_omega/(4*pi)));        % coarse mesh size


% width of PML
high_wpml = 8*high_wl(1)*ones(size(high_omega)); %fh.*ceil(high_wl./fh);
low_wpml = 2*ch.*ceil(low_wl(1)./ch);


%% Generate the domain sizes
sd = 1/2;
Rest = 1.5;           % estimate of the distance to the source point

high_r = NMLA_radius(high_omega,Rest);
md = sd + high_r + high_wpml;
md = ceil(md*10)/10;      % middle domain size

% Rest = sqrt(2)*md;
low_r = NMLA_radius(low_omega,Rest);
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
    A = assemble_Helmholtz_matrix_SFEM(lnode,lelem,omega,wpml,sigmaMax,speed,fquadorder);
    b = assemble_RHS_SFEM(lnode,lelem, @(x)nodal_basis(xs,ys,h,x),fquadorder);
    b = b/(h*h/2);
    [~,~,isBdNode] = findboundary(lelem);
    freeNode = find(~isBdNode);
    lN = size(lnode,1);        u_std = zeros(lN,1);
    u_std(freeNode) = A(freeNode,freeNode)\b(freeNode);
    toc;
    
    
    %% Step 2: Use NMLA to find ray directions d_c with low frequency sqrt(\omega)
    fprintf(['-'*ones(1,80) '\n']);
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
        if r0>2*epsilon
            [cnumray_angle(i)] = NMLA(x0,y0,c0,omega,Rest,lnode,lelem,u_std,ux,uy,[],1/5,Nray,'num',sec_opt,plt);
        else
            cnumray_angle(i) =  ex_ray([x0,y0],xs,ys,0);
        end
    end
    cnumray = exp(1i*cnumray_angle);
    ray = interpolation(cnode,celem,mnode,cnumray);
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
        if r0>2*epsilon
            %             Rest = r0;
            [cnumray_angle(i)] = NMLA(x0,y0,c0,omega,Rest,mnode,melem,uh1,ux,uy,[],1/5,Nray,'num',sec_opt,plt);
        else
            cnumray_angle(i) =  ex_ray([x0,y0],xs,ys,0);
        end
    end
    cnumray = exp(1i*cnumray_angle);
    ray = interpolation(cnode,celem,node,cnumray);
    toc;
    
    
    %% Step 5: Solve the original Helmholtz equation by Ray-based FEM with ray directions d_o
    fprintf([ '-'*ones(1,80) '\n']);
    fprintf('Step5: Ray-FEM, high frequency \n');
    tic;
    
    % Assembling
    omega = high_omega(ti);
    wpml = high_wpml(ti);                % width of PML
    sigmaMax = 25/wpml;                 % Maximun absorbtion
    
    option = 'homogeneous';
    A = assemble_Helmholtz_matrix_RayFEM(node,elem,omega,wpml,sigmaMax,speed,ray,fquadorder);
    b = assemble_RHS_RayFEM_with_ST(node,elem,xs,ys,omega,epsilon,wpml,sigmaMax,ray,speed,fquadorder,option);
    [~,v] = RayFEM_direct_solver(node,elem,A,b,omega,ray,speed);
    toc;
    
    clear lnode lelem mnode melem cnode celem cnumray cnumray_angle;
    clear A b x y rr cf u_std ub uh uh1 ux uy freeNode isBdNode;
    
    
    %% Compute errors
%     fprintf([ '-'*ones(1,80) '\n']);
%     fprintf('Compute errors\n');
%     tic;
%     switch round(omega/(2*pi))
%         case 60
%             load('../Solutions_Lippmann_Schwinger/point_source_k_60_2pi.mat')
%         case 80
%             load('../Solutions_Lippmann_Schwinger/point_source_k_80_2pi.mat')
%         case 120
%             load('../Solutions_Lippmann_Schwinger/point_source_k_120_2pi.mat')
%         case 160
%             load('../Solutions_Lippmann_Schwinger/point_source_k_160_2pi.mat')
%         case 240
%             load('../Solutions_Lippmann_Schwinger/point_source_k_240_2pi.mat')
%     end
% 
%     rh = 1/4000;
%     [rnode,~] = squaremesh([-a,a,-a,a],rh); rn = round(sqrt(size(rnode,1)));
%     uh = RayFEM_solution(node,elem,omega,speed,v,ray,rnode);
%     
%     % Reference solution 
%     x = rnode(:,1); y = rnode(:,2);
%     rr = sqrt((x-xs).^2 + (y-ys).^2);
%     
%     cf = cutoff(epsilon,2*epsilon,rnode,xs,ys);
%     ur = (1-cf).*u;
%     ur(rr<epsilon) = 0;
%    
%     % Errors
%     du = uh - ur;
%     idx = find( ~( (x<=max(x)-wpml).*(x>= min(x)+wpml)...
%         .*(y<= max(y)-wpml).*(y>= min(y)+wpml) ) ); % index on PML
%     du(idx) = 0;  ur(idx) = 0;
%     
%     max_err(ti) = norm(du,inf);
%     rel_max_err(ti) = norm(du,inf)/norm(ur,inf);
%     l2_err(ti) = norm(du)*h;
%     rel_l2_err(ti) = norm(du)/norm(ur);
%     toc;
%     
%     clear node elem rnode x y rr cf idx u;
    
    
    %% Plots
%     sh = 1/400;
%     [snode,selem] = squaremesh([-a,a,-a,a],sh); 
%     sn = round(sqrt(size(snode,1))); idx = 1:sn; idx = 10*(idx-1) + 1;
%     su = reshape(du,rn,rn); su = su(idx,idx);
%     figure(8); 
%     subplot(1,2,1);
%     showsolution(snode,selem,real(su(:))); colorbar;
%     title('Ray-FEM solution error')
%     subplot(1,2,2);
%     showsolution(snode,selem,real(su(:)),2); colorbar;
%     title('Ray-FEM solution error')
%     
%     figure(9);
%     subplot(1,2,1);
%     showsolution(snode,selem,speed(snode)); colorbar;
%     title('Wave speed')
%     subplot(1,2,2);
%     showsolution(snode,selem,speed(snode),2); colorbar;
%     title('Wave speed')
    
end

totaltime = toc(tstart);
fprintf('\n\nTotal running time: % d minutes \n', totaltime/60);



%% plots
% figure(2);
% subplot(2,2,1);
% show_convergence_rate(high_omega(1:test_num),max_err(1:test_num),'omega','max err');
% subplot(2,2,2);
% show_convergence_rate(high_omega(1:test_num),l2_err(1:test_num),'omega','L2 err');
% subplot(2,2,3);
% show_convergence_rate(high_omega(1:test_num),rel_max_err(1:test_num),'omega','Rel max ');
% subplot(2,2,4);
% show_convergence_rate(high_omega(1:test_num),rel_l2_err(1:test_num),'omega','Rel L2 ');

figure(3);
subplot(1,2,1);
show_convergence_rate(high_omega(1:test_num),l2_err(1:test_num),'omega','L2 err');
subplot(1,2,2);
show_convergence_rate(high_omega(1:test_num),rel_l2_err(1:test_num),'omega','Rel L2 err');



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
