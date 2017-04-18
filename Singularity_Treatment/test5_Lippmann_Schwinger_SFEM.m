%% Test the convergence order of SFEM for singularity removel problem: error \sim O(h^2)

% Set up
omega = 120*pi;
wl = 2*pi/omega;
epsilon = 1/(2*pi);

%% Load source/wavespeed data
xs = -0.2; ys = -0.2;                     % point source location

sigma = 0.15;
xHet = 0.2;   yHet = 0.2;

nu = @(x,y) 0.2*exp( -1/(2*sigma^2)*((x-xHet).^2 + (y-yHet).^2) )...
    .*Lippmann_Schwinger_window(sqrt((x-xHet).^2 + (y-yHet).^2), 0.22,0.16  );

speed = @(p) 1./sqrt(1 + nu( p(:,1), p(:,2) ));    % wave speed
speed_min = 1/sqrt(1.2);


wpml = 0.1;
sigmaMax = 25/wpml;
fquadorder = 3;
a = 1/2;

nt = 3;
% hs = zeros(nt,1);
SFEM_errs = zeros(nt,1);  % error on physical domain
hs = 1./[1800 2000 2400]';
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
    [bdNode,~,isBdNode] = findboundary(elem);
    freeNode = find(~isBdNode);
    
    % S-FEM solution
    N = size(node,1);  n = round(sqrt(N));
    us = zeros(N,1);
    us(freeNode) = A(freeNode,freeNode)\b(freeNode);
    
    % Exact solution
    load('/Solutions_Lippmann_Schwinger/point_source_k_120_jun_new_version.mat');
    rh = 1/5000;
    [rnode,relem] = squaremesh([-a,a,-a,a],rh); rn = round(sqrt(size(rnode,1)));
    ur = interpolation(rnode,relem,node,u);
    
    x = node(:,1); y = node(:,2);
    rr = sqrt((x-xs).^2 + (y-ys).^2);
%     ub = 1i/4*besselh(0,1,omega*rr);
    cf = cutoff(epsilon,2*epsilon,node,xs,ys);
    uex = (1-cf).*ur;
    uex(rr<epsilon) = 0;
    
    % Error
    du = us - uex;
    idx = find( (x<=max(x)-wpml).*(x>= min(x)+wpml)...
        .*(y<= max(y)-wpml).*(y>= min(y)+wpml) );
    du_phy = du(idx);
    SFEM_errs(ii) = norm(du_phy)/norm(uex(idx));
    toc;
end

% Plot
figure(11);
show_convergence_rate(hs,SFEM_errs,'h','||u - u_h||_{L^2(\Omega)}');





