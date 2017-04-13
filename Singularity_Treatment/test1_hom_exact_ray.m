%% Convergence test for homogenenous case with exact ray information

clear;
addpath(genpath('../../ifem/'));
addpath('../Methods/');
addpath('../Cutoff_Functions/')
addpath('../Plots_Prints/');


%% Set up
% xs = -0.2; ys = -0.2;                     % source location 
xs = 0; ys = 0;                     % source location 
epsilon = 1/(2*pi);               % cut-off parameter   
speed = @(p) ones(size(p(:,1)));    % wave speed

wpml = 0.1;                         % width of PML  
sigmaMax = 25/wpml;                 % absorption of PML  
fquadorder = 3;                     % numerical quadrature order 
a = 1/2;                            % computational domain [-a,a]^2

NPW = 4;                            % grid number per wavelength
omegas = pi*[120 160 240 320 480 640];      % omega's
nt = 5; length(omegas);                % number of tests
rel_l2_err = zeros(1,nt);           % record relative l2 errors


%% Tests
for ni = 1:nt
    
    omega = omegas(ni);
    wl = 2*pi/omega;
    h = 1/round(1/(wl/NPW));
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

%% save output 
nameFile = strcat('resutls_1_HomExRay_NPW_', num2str(NPW), '.mat');
save(nameFile, 'rel_l2_err', 'NPW', 'omegas', 'test_num');


%% plot
figure(22);
show_convergence_rate(omegas(1:nt), rel_l2_err,'omega','||u - u_h||_{L^2(\Omega)}/||u||_{L^2(\Omega)}');



