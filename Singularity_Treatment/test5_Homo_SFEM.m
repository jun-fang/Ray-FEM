%% Test the convergence order of SFEM for singularity removel problem: error \sim O(h^2)

clear;
% Set up
omega = 120*pi;
wl = 2*pi/omega;
epsilon = 1/(2*pi);

%% Load source/wavespeed data
xs = -0.2; ys = -0.2;                     % point source location

% sigma = 0.15;
% xHet = 0.2;   yHet = 0.2;
% 
% nu = @(x,y) 0.2*exp( -1/(2*sigma^2)*((x-xHet).^2 + (y-yHet).^2) )...
%     .*Lippmann_Schwinger_window(sqrt((x-xHet).^2 + (y-yHet).^2), 0.22,0.16  );

% speed = @(p) 1./sqrt(1 + nu( p(:,1), p(:,2) ));    % wave speed
% speed_min = 1/sqrt(1.2);

speed = @(p) 1./ones(size(p(:,1)));

wpml = 0.1;
sigmaMax = 25/wpml;
fquadorder = 3;
a = 1/2;

nt = 5;
% hs = zeros(nt,1);
SFEM_errs = zeros(nt,1);  % error on physical domain
hs = 1./[1200 1500 1800 2100 2400]';
hs = hs(1:nt);

% Tests
for ii = 1:nt
    tic;
    h = hs(ii);
%     h = h/2;
%     hs(ii) = h;
    [node,elem] = squaremesh([-a,a,-a,a],h);
    fprintf('Case %d: omega/pi = %d, 1/h = %d\n', ii, round(omega/pi), round(1/h));

    % Aseemble
    A = assemble_Helmholtz_matrix_SFEM(node,elem,omega,wpml,sigmaMax,speed,fquadorder);
    b = assemble_RHS_SFEM_with_ST(node,elem,xs,ys,omega,wpml,sigmaMax,epsilon,fquadorder);
    
    % Boundary conditions
    [~,~,isBdNode] = findboundary(elem);
    freeNode = find(~isBdNode);
    
    % S-FEM solution
    N = size(node,1);  n = round(sqrt(N));
    us = zeros(N,1);
    us(freeNode) = A(freeNode,freeNode)\b(freeNode);
    
%     % Exact solution
%     load('/Solutions_Lippmann_Schwinger/point_source_k_120_jun_new_version.mat');
%     rh = 1/5000;
%     [rnode,relem] = squaremesh([-a,a,-a,a],rh); rn = round(sqrt(size(rnode,1)));
%     ur = interpolation(rnode,relem,node,u);
    
    x = node(:,1); y = node(:,2);
    rr = sqrt((x-xs).^2 + (y-ys).^2);
    ub = 1i/4*besselh(0,1,omega*rr);
    cf = cutoff(epsilon,2*epsilon,node,xs,ys);
    uex = (1-cf).*ub;
    uex(rr<epsilon) = 0;
    
    % Error
    du = us - uex;
    idx = find( ~((x<=max(x)-wpml).*(x>= min(x)+wpml)...
        .*(y<= max(y)-wpml).*(y>= min(y)+wpml) ) );
    du(idx) = 0;
    SFEM_errs(ii) = norm(du)/norm(uex(idx));
    toc;
end

clear A b x y rr ub cf freeNode idx isBdNode

% Plot
figure(11);
show_convergence_rate(hs,SFEM_errs,'h','||u - u_h||_{L^2(\Omega)}');

% print results
fprintf( ['\n' '-'*ones(1,80) '\n']);
fprintf( 'NPW:                 ');
fprintf( '&  %.2e  ',round(wl./hs) );
fprintf( '\n1/h:                 ');
fprintf( '&  %.2e  ',round(1./hs));
fprintf( '\nSFEM rel l2 err:     ');
fprintf( '&  %.2e  ',SFEM_errs );
fprintf( ['\n' '-'*ones(1,80) '\n']);


% --------------------------------------------------------------------------------
% NPW:                 &  2.00e+01  &  2.50e+01  &  3.00e+01  &  3.50e+01  &  4.00e+01  
% 1/h:                 &  1.20e+03  &  1.50e+03  &  1.80e+03  &  2.10e+03  &  2.40e+03  
% SFEM rel l2 err:     &  8.02e-01  &  5.43e-01  &  3.85e-01  &  2.86e-01  &  2.20e-01  
% --------------------------------------------------------------------------------
figure(12);
subplot(1,3,1);
showsolution(node,elem,real(uex),2);colorbar;
title('Homogeneous exact solu: 120pi');

subplot(1,3,2);
showsolution(node,elem,real(us),2);colorbar;
title('Homo SFEM solu: 120pi, NPW = 40');

subplot(1,3,3);
showsolution(node,elem,real(du),2);colorbar;
title('Homo SFEM error: 120pi, NPW = 40');


