%% Convergence test for gravity Helmholtz equation with exact rays

clear;

addpath(genpath('../../ifem/'));
addpath(genpath('../../lhelmfs/'));
addpath('../Methods/');
addpath('../NMLA/');
addpath('../Cutoff_Functions/')
addpath('../Plots_Prints/');


xs = 0; ys = 0;                     % source location
epsilon = 50/(80*pi);               % cut-off parameter
omega0 = 400*pi;
E0 = omega0^2;                   % E0 should be large enough to make sure the medium is smooth enough
speed = @(p) omega0./sqrt( E0 + p(:,2) );    % wave speed


node0 = [xs, ys];
S02 = ( 1/speed(node0) )^2;
g0 = [0, 1/(2*omega0*omega0)];


wpml = 0.1;                         % width of PML
sigmaMax = 25/wpml;                 % absorption of PML
fquadorder = 3;                     % numerical quadrature order
a = 1/2;                            % computational domain [-a,a]^2


nt = 4;                             % number of tests
errors = zeros(1,nt);
rhss = zeros(1,nt);
omegas = pi*[120,160,240,320];      % omega's
NPW = 4;                            % grid number per wavelength


for ii = 1:nt
    tic;
    
    % omega and mesh
    omega = omegas(ii);
    wl = 2*pi/omega;
    h = 1/round(1/(wl/NPW));
    [node,elem] = squaremesh([-a,a,-a,a],h);
    fprintf('\nCase %d: omega/pi = %d,  NPW = %d,  1/h = %d\n', ii, round(omega/pi), NPW, 1/h);
    
    % Exact ray information
    [dx, dy] = eikonal_cgss(S02, g0, node0, node);
    dr2 = dx.^2 + dy.^2;
    ray_angle = atan2(dy,dx);
    ray = exp(1i*ray_angle).*(dr2>10*eps);
    
    % Gravity parameters
    alpha = (omega/omega0)^2;
    E = E0*alpha;
    
    % Assemble and solve the system Au = b
    option = {'gravity', alpha, E};
    A = assemble_Helmholtz_matrix_RayFEM(node,elem,omega,wpml,sigmaMax,speed,ray,fquadorder);
    b = assemble_RHS_RayFEM_with_ST(node,elem,xs,ys,omega,epsilon,wpml,sigmaMax,ray,speed,fquadorder,option);
    u = singularity_direct_solver(node,elem,A,b,omega,ray,speed);
    
    
    % get the exact solution
    trg = node';  src = [xs;ys];
    ub = lhelmfs(trg,src,alpha,E,0);
    ub = ub(:);
    
    cf = cutoff(epsilon,2*epsilon,node,xs,ys);
    uex = (1-cf).*ub;
    
    x = node(:,1);  y = node(:,2);
    xx = x - xs; yy = y - ys;
    rr = sqrt(xx.*xx + yy.*yy);
    uex(rr<epsilon) = 0;
    
    % compute the error
    du = u - uex;
    idx = find( ~( (x<=max(x)-wpml).*(x>= min(x)+wpml)...
        .*(y<= max(y)-wpml).*(y>= min(y)+wpml) ) ); % index on PML
    du(idx) = 0;  uex(idx) = 0;
    errors(ii) = norm(du)/norm(uex);
    
    toc;
end


%% plot
figure(33);
hold off
show_convergence_rate(omegas(1:nt),errors,'omega','||u - u_h||_{L^2(\Omega)}/||u||_{L^2(\Omega)}');


