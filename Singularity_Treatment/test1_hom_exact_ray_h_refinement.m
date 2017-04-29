%% h convergence test for homogenenous case with exact ray information

clear;
addpath(genpath('../../ifem/'));
addpath('../Methods/');
addpath('../Functions/')
addpath('../Plots_Prints/');


%% Set up
xs = 0; ys = 0;                     % source location 
epsilon = 1/(2*pi);                 % cut-off parameter   
speed = @(p) ones(size(p(:,1)));    % wave speed

wpml = 0.1;                         % width of PML  
sigmaMax = 25/wpml;                 % absorption of PML  
fquadorder = 3;                     % numerical quadrature order 
a = 1/2;                            % computational domain [-a,a]^2

NPW = 4;                            % grid number per wavelength
test_num = 3;                       % number of tests
rel_l2_err = zeros(1,test_num); 

% frequency and mesh size
omega = 100*pi;
wl = 2*pi/omega;
h = 2/round(1/(wl/NPW));
hs = h./[2 4 8];


%% Tests
for ni = 1:test_num
    
    h = hs(ni);
    [node,elem] = squaremesh([-a,a,-a,a],h);
    fprintf(['-'*ones(1,80) '\n']);
    fprintf('Case %d: omega/pi = %d, NPW = %d, 1/h = %d\n', ni, round(omega/pi), NPW, round(1/h));
    
    tic;
    %% Exact ray information
    x = node(:,1);  y = node(:,2);
    rr = sqrt((x-xs).^2 + (y-ys).^2);
    ray = atan2(y-ys, x-xs);
    ray = exp(1i*ray).*(rr>10*eps);
    
    %% Assemble and solve the system Au = b
    option = 'homogeneous';
    A = assemble_Helmholtz_matrix_RayFEM(node,elem,omega,wpml,sigmaMax,speed,ray,fquadorder);
    b = assemble_RHS_RayFEM_with_ST(node,elem,xs,ys,omega,epsilon,wpml,sigmaMax,ray,speed,fquadorder,option);
    u = RayFEM_direct_solver(node,elem,A,b,omega,ray,speed);
       
    %% Get the exact solution
    ub = 1i/4*besselh(0,1,omega*rr);
    cf = cutoff(epsilon,2*epsilon,node,xs,ys);
    uex = (1-cf).*ub;
    uex(rr<epsilon) = 0;
    
    %% Compute relative L2 error 
    du = u - uex;
    idx = find( ~( (x<=max(x)-wpml).*(x>= min(x)+wpml)...
        .*(y<= max(y)-wpml).*(y>= min(y)+wpml) ) ); % index on PML
    du(idx) = 0;  uex(idx) = 0;
    rel_l2_err(ni) = norm(du)/norm(uex);
    toc;
end


%% plot
figure(12);
show_convergence_rate(hs(1:test_num), rel_l2_err,'h','||u - u_h||_{L^2(\Omega)}/||u||_{L^2(\Omega)}');



