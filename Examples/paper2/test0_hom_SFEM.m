%% Test the convergence order of SFEM for singularity removel problem: error \sim O(h^2)

%% Set path
addpath(genpath('../../../ifem/'));
addpath('../../Methods/');
addpath('../../Functions/')
addpath('../../Plots_Prints/');
addpath('../../Singularity_Treatment/');


%% Set up
xs = 0; ys = 0;
omega = 40*pi;
wl = 2*pi/omega;
epsilon = 0.143;
speed = @(p) ones(size(p(:,1)));

wpml = 0.2;
sigmaMax = 25/wpml;
fquadorder = 3;
a = 1/2;

nt = 3;
hs = zeros(nt,1);
SFEM_errs = zeros(nt,1);  % error on physical domain
h = 1/100;

%% Tests
for ii = 1:nt
    tic;
    h = h/2;
    hs(ii) = h;
    [node,elem] = squaremesh([-a,a,-a,a],h);
    fprintf('Case %d: omega/pi = %d, 1/h = %d\n', ii, round(omega/pi), round(1/h));

    % Aseemble
    A = assemble_Helmholtz_matrix_SFEM(node,elem,omega,wpml,sigmaMax,speed,fquadorder);
    b = assemble_RHS_SFEM_with_ST(node,elem,xs,ys,omega,wpml,sigmaMax,epsilon,fquadorder);
    
    % Boundary conditions
    [bdNode,~,isBdNode] = findboundary(elem);
    freeNode = find(~isBdNode);
    
    % S-FEM solution
    N = size(node,1);  n = round(sqrt(N));
    u = zeros(N,1);
    u(freeNode) = A(freeNode,freeNode)\b(freeNode);
    
    % Exact solution
    x = node(:,1); y = node(:,2);
    rr = sqrt((x-xs).^2 + (y-ys).^2);
    ub = 1i/4*besselh(0,1,omega*rr);
    cf = cutoff(epsilon,2*epsilon,node,xs,ys);
    uex = (1-cf).*ub;
    uex(rr<epsilon) = 0;
    
    % Error
    du = u - uex;
    idx = find( (x<=max(x)-wpml).*(x>= min(x)+wpml)...
        .*(y<= max(y)-wpml).*(y>= min(y)+wpml) );
    du_phy = du(idx);
    SFEM_errs(ii) = norm(du_phy)/norm(uex(idx));
    toc;
end

% Plot
figure(11);
show_convergence_rate(hs,SFEM_errs,'h','||u - u_h||_{L^2(\Omega)}');





