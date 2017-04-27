%% Test the convergence order of SFEM for singularity removel problem: error \sim O(h^2)


clear;

addpath(genpath('../../ifem/'));
addpath(genpath('../../lhelmfs/'));
addpath('../Methods/');
addpath('../NMLA/');
addpath('../Cutoff_Functions/')
addpath('../Plots_Prints/');
    

xs = 0; ys = 0;                     % source location
epsilon = 0.379; 1/(2*pi);               % cut-off parameter
omega0 = 10; 0.8;
E0 = omega0^2;                   % E0 should be large enough to make sure the medium is smooth enough
speed = @(p) omega0./sqrt( E0 + p(:,2) );    % wave speed


node0 = [xs, ys];
S02 = ( 1/speed(node0) )^2;
g0 = [0, 1/(2*omega0*omega0)];


omega = 40*pi;
wl = 2*pi/omega;

wpml = 0.1;
sigmaMax = 25/wpml;
fquadorder = 3;
a = 0.9;

nt = 4;
% hs = zeros(nt,1);
SFEM_errs = zeros(nt,1);  % error on physical domain
hs = 1./[400 600 800 1000 ]';
% hs = 1./[800  1000  1200  1600  1800]';
hs = hs(1:nt);

% Tests
for ii = 1:nt
    tic;
    h = hs(ii);
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
    
    % get the exact solution
    alpha = (omega/omega0)^2;
    E = alpha*E0;
    
    trg = node';  src = [xs;ys];
    ub = lhelmfs(trg,src,alpha,E,0);
    ub = ub(:);
    
    cf = cutoff(epsilon,2*epsilon,node,xs,ys);
    uex = (1-cf).*ub;
    
    x = node(:,1);  y = node(:,2);
    xx = x - xs; yy = y - ys;
    rr = sqrt(xx.*xx + yy.*yy);
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
% figure(12);
% subplot(1,3,1);
% showsolution(node,elem,real(uex),2);colorbar;
% title('Homogeneous exact solu: 120pi');
% 
% subplot(1,3,2);
% showsolution(node,elem,real(us),2);colorbar;
% title('Homo SFEM solu: 120pi, NPW = 40');
% 
% subplot(1,3,3);
% showsolution(node,elem,real(du),2);colorbar;
% title('Homo SFEM error: 120pi, NPW = 40');
% 
% 
