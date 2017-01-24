%% Convergence test for gravity Helmholtz equation with exact ray information

clear;

addpath(genpath('../../ifem/'));
addpath(genpath('../../lhelmfs/'));
addpath('../Methods/');
addpath('../NMLA/');
addpath('../Cutoff_Functions/')
addpath('../Plots_Prints/');


xs = 0; ys = 0;                     % source location
epsilon = 50/(80*pi);               % cut-off parameter
omega0 = 80*pi;
E0 = omega0^2;
speed = @(p) omega0./sqrt(E0 + (p(:,1)));    % wave speed

S02 = E0/omega0^2;
g0 = [0, 1/(2*omega0*omega0)];
node0 = [xs, ys];


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

    omega = omegas(ii);
    wl = 2*pi/omega;
    h = 1/round(1/(wl/NPW));
    [node,elem] = squaremesh([-a,a,-a,a],h);
    fprintf('Case %d: omega/pi = %d,  1/h = %d\n', ii, round(omega/pi), 1/h);
    
    % Exact ray information
    [dx, dy] = eikonal_cgss(S02, g0, node0, node);
    dr2 = dx.^2 + dy.^2;
    ray_angle = atan2(dy,dx);
    ray = exp(1i*ray_angle).*(dr2>10*eps);
    
    % Gravity parameters
    alpha = (omega/omega0)^2;
    E = E0*alpha;
    
    % Ray-FEM assembling with singularity treatment
    tic;
    A = assemble_Helmholtz_matrix_with_ray_1(node,elem,omega,wpml,sigmaMax,speed,ray,fquadorder);
    b = assemble_RHS_gravity(node,elem,xs,ys,epsilon,wpml,sigmaMax,omega,E,alpha,speed,ray,fquadorder);
    
    % Boundary conditions
    [~,~,isBdNode] = findboundary(elem);
    freeNode = find(~isBdNode);
    N = size(node,1);  Nray = size(ray,2);
    v = zeros(N,1);
    v(freeNode) = A(freeNode,freeNode)\b(freeNode);
    
    % construct solution
    grad = ray(:);
    grad = [real(grad),imag(grad)];
    repnode = repmat(node,Nray,1);
    temp = grad(:,1).*repnode(:,1) + grad(:,2).*repnode(:,2);
    
    k = omega./speed(node);           % wavenumber
    kk = repmat(k,1,Nray);
    u = v.*exp(1i*kk(:).*temp);
    u = reshape(u,N,Nray);
    u = sum(u,2);
    
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
    du = u - uex;
    
    idx = find( (x<=max(x)-wpml).*(x>= min(x)+wpml)...
        .*(y<= max(y)-wpml).*(y>= min(y)+wpml) );
    
    dd = 0*du; dd(idx) = du(idx);
    
    errors(ii) = norm(dd)*h; %/norm(uex(idx));
    toc;
end


%% plot
figure(33);
show_convergence_rate(omegas(1:nt),errors,'omega','||u - u_h||_{L^2(\Omega)}');



