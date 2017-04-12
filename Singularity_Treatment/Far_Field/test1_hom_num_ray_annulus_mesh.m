%% Compute errors of Ray-FEM solution in homogeneous medium
function [] = test1_hom_num_ray_annulus_mesh(NPW, test_num)

addpath(genpath('../../../ifem/'));
addpath('../../Methods/');
addpath('../../NMLA/');
addpath('../../Cutoff_Functions/')
addpath('../../Plots_Prints/');
addpath('../../Singularity_Treatment/');

xs = 0;   ys = 0;                  % point source location
speed = @(x) ones(size(x,1),1);    % medium speed

plt = 0;                   % show solution or not
sec_opt = 0;               % NMLA second order correction or not
Rest = 0.35;

fquadorder = 3;                    % numerical quadrature order
epsilon = 50/(80*pi);
% NPW = 8;
wpml = 0.075;%3*round(wavelength/h)*h;                 % width of PML
sigmaMax = 25/wpml;                % Maximun absorbtion

a = 0.5; bd = 0.1;

% test_num = 3;
omegas = [120 160 200 240 320 480 640 800 960]*pi;
high_omega = omegas;

fh = 1./(NPW*round(omegas/(2*pi)));      % fine mesh size
% ch = 1/10*ones(size(fh));

ch = 1./(10*round(sqrt(2./fh)/10));
% ch = 1./(20*round(sqrt(omegas)/(2*pi)));

% ch = fh; %1/20*ones(size(fh));

% error
high_max_rayerr = 0*high_omega;    % L_inf ray error of high-freq waves
high_l2_rayerr = 0*high_omega;     % L_2 ray error of high-freq waves

max_err = 0*high_omega;            % L_inf error of the numerical solution
rel_max_err = 0*high_omega;        % relative L_inf error of the numerical solution
l2_err = 0*high_omega;             % L_2 error of the numerical solution
rel_l2_err = 0*high_omega;         % relative L_2 error of the numerical solution


tstart = tic;
for ti = 1:test_num
    
    omega = omegas(ti);
    h = fh(ti);  h_c = ch(ti);
    fprintf(['-'*ones(1,80) '\n']);
    fprintf('\nCase %d: \n  omega/(2*pi) = %d,   1/h = %d,  1/h_c = %d,  NPW = %d \n',...
        ti, round(omega/(2*pi)), 1/h, 1/h_c, NPW);
    
    
    [node,elem] = squaremesh_annulus([-a,a,-a,a],[-bd,bd,-bd,bd],h);
    Nray = 1; N = size(node,1); Ndof = N;
    
    [lnode,lelem] = squaremesh([-0.6,0.6,-0.6,0.6],h);
    x = lnode(:,1); y = lnode(:,2);
    rr = sqrt((x-xs).^2 + (y-ys).^2);
    ub = 1i/4*besselh(0,1,omega*rr); 
    [ux,uy] = num_derivative(ub,h,2);
    
    %% Exact ray information
    x = node(:,1);  y = node(:,2);
    rr = sqrt((x-xs).^2 + (y-ys).^2);
    ray = atan2(y-ys, x-xs);
    exray = exp(1i*ray).*(rr>10*eps);
    
    
    %% NMLA
    fprintf([ '-'*ones(1,80) '\n']);
    fprintf('Step2: NMLA, low frequency \n');
    
    [cnode,celem] = squaremesh([-a,a,-a,a],h_c);
    cN = size(cnode,1);
    cnumray_angle = zeros(cN,Nray);
    
    % NMLA
    tic;
    for i = 1:cN
        x0 = cnode(i,1);  y0 = cnode(i,2);
        r0 = sqrt((x0-xs)^2 + (y0-ys)^2);
        c0 = speed(cnode(i,:));
        if r0 > bd - 3*h_c
%             Rest = r0;
            [cnumray_angle(i)] = NMLA(x0,y0,c0,omega,Rest,lnode,lelem,ub,ux,uy,[],1/5,Nray,'num',sec_opt,plt);
        end
    end
    cnumray = exp(1i*cnumray_angle);
    numray = interpolation(cnode,celem,node,cnumray);
    numray = numray./abs(numray);
    numray = numray.*(rr > 2*epsilon) + exray.*(rr <= 2*epsilon);
    toc;
    
    rayerr = numray - exray;
    high_max_rayerr(ti) = norm(rayerr,inf);
    high_l2_rayerr(ti) = norm(rayerr)*h/(norm(exray)*h);    
    
    
    
    %% Assemble and solve the system Au = b
    fprintf(['-'*ones(1,80) '\n']);
    fprintf('Step3: Ray-FEM, high frequency \n');
    
    tic;
    
    ray = numray;
    option = 'homogeneous';
    A = assemble_Helmholtz_matrix_RayFEM(node,elem,omega,wpml,sigmaMax,speed,ray,fquadorder);
    b = assemble_RHS_RayFEM_with_ST(node,elem,xs,ys,omega,epsilon,wpml,sigmaMax,ray,speed,fquadorder,option);
    u = RayFEM_direct_solver(node,elem,A,b,omega,ray,speed);
    
    
    %% Get the exact solution
    ub = 1i/4*besselh(0,1,omega*rr);
    cf = cutoff(epsilon,2*epsilon,node,xs,ys);
    uex = (1-cf).*ub;
    uex(rr<epsilon) = 0;
    
    
    %% Errors
    du = u - uex;
    idx = find( ~( (x <= max(x)-wpml).*(x >= min(x)+wpml)...
        .*(y <= max(y)-wpml).*(y >= min(y)+wpml) ) ); % index on PML
    du(idx) = 0;  uex(idx) = 0;
    
    max_err(ti) = norm(du,inf);
    rel_max_err(ti) = norm(du,inf)/norm(uex,inf);
    l2_err(ti) = norm(du)*h;
    rel_l2_err(ti) = norm(du)/norm(uex);
    
    toc;
    
end


totaltime = toc(tstart);
fprintf('\n\nTotal running time: % d minutes \n', totaltime/60);


%% save output 
nameFile = strcat('resutls_1_HomNumRay_NPW_', num2str(NPW), '.mat');
save(nameFile, 'rel_l2_err', 'NPW', 'high_omega');




%% plots
figure(1);
subplot(1,2,1);
show_convergence_rate(high_omega(1:test_num),high_max_rayerr(1:test_num),'omega','high max');
subplot(1,2,2);
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
