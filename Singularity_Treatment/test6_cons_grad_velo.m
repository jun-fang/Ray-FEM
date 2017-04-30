%% h convergence test for homogenenous case with exact ray information

clear;
addpath(genpath('../../ifem/'));
addpath('../Methods/');
addpath('../Functions/')
addpath('../Plots_Prints/');


%% Set up
xs = 0; ys = 0;                     % source location 
epsilon = 1/(2*pi);                 % cut-off parameter   

v0 = 1; G0 = [0.1, -0.2];
speed = @(p) 1./( v0 + (G0(1)*p(:,1) + G0(2)*p(:,2)) );    % wave speed


wpml = 0.1;                         % width of PML  
sigmaMax = 25/wpml;                 % absorption of PML  
fquadorder = 3;                     % numerical quadrature order 
a = 1/2;                            % computational domain [-a,a]^2

NPW = 4;                            % grid number per wavelength
test_num = 7;                       % number of tests
rel_l2_err = zeros(1,test_num); 

% frequency and mesh size
omega = 100*pi;
wl = 2*pi/omega;
% h = 2/round(1/(wl/NPW));
hs = 1./[200 280 350 400 560 700 1000];
rn = [10 7 5 4];
rh = 1/2800;
[rnode, ~] = squaremesh([-a, a, -a, a], rh);
x = rnode(:,1);  y = rnode(:,2);
rr = sqrt((rnode(:,1)-xs).^2 + (rnode(:,2)-ys).^2);
load('/Users/junfang/Dropbox/test6_reference_solution.mat');


%% Babich expansion

load('Babich_CGV.mat');

CompressRatio = 20;
Bh = 1/round( 1/(Bh0/CompressRatio) );
Bx = -a: Bh : a;  By = -a: Bh : a;
[BX0, BY0] = meshgrid(Bx0, By0);
[BX, BY] = meshgrid(Bx, By);


%% refined phase and amplitude
% ttao = interp2(BX0,BY0,tao,BX,BY,'spline'); % refined phase
DD1 = interp2(BX0,BY0,D1,BX,BY,'spline'); % refined amplitude
DD2 = interp2(BX0,BY0,D2,BX,BY,'spline'); % refined amplitude


%% analytical ray and phase
[ray, T, ~, ~] = eikonal_cgv(1, [0.1, -0.2], [0,0], [BX(:), BY(:)]);
T = reshape(T, size(BX));

%% Babich expansion
G1 = 1i*(sqrt(pi)/2)*besselh(0,1,omega*T);
G1 = G1.*DD1;

G2 = 1i*(sqrt(pi)/2)*exp(1i*pi)*(2*T/omega);
G2 = G2.*besselh(1,1,omega*T);
G2 = G2.*DD2;

ub = G1 + G2;
ub = ub(:);

clear BX BY BX0 BY0 Bx By Bx0 By0 D1 D2 DD1 DD2 G1 G2;


%% Tests
for ni = 1:test_num
    
    h = hs(ni);
    [node,elem] = squaremesh([-a,a,-a,a],h);
    fprintf(['-'*ones(1,80) '\n']);
    fprintf('Case %d: omega/pi = %d, NPW = %d, 1/h = %d\n', ni, round(omega/pi), NPW, round(1/h));
    
    tic;
    option ={ 'Babich_CGV', 'exact_phase'};
    ray = eikonal_cgv(1, [0.1, -0.2], [0,0], node);
    A = assemble_Helmholtz_matrix_RayFEM(node,elem,omega,wpml,sigmaMax,speed,ray,fquadorder);
    b = assemble_RHS_RayFEM_with_ST(node,elem,xs,ys,omega,epsilon,wpml,sigmaMax,ray,speed,fquadorder,option);
    [~, v] = RayFEM_direct_solver(node,elem,A,b,omega,ray,speed);
    uh = RayFEM_solution(node,elem,omega,speed,v,ray,rnode);
    toc;
    
    %% singularity part
    cf = cutoff(epsilon,2*epsilon,rnode,xs,ys);
    u = uh + cf.*ub;

    
    %% Compute relative L2 error 
    du = u - u_ref;
    idx = find( ~( (x<=max(x)-wpml).*(x>= min(x)+wpml)...
        .*(y<= max(y)-wpml).*(y>= min(y)+wpml).*(rr>0.1) ) ); % index on PML
    du(idx) = 0;  u_ref(idx) = 0;
    rel_l2_err(ni) = norm(du)/norm(u_ref);
    toc;
end


%% plot
figure(12);
show_convergence_rate(hs(1:test_num), rel_l2_err,'h','||u - u_h||_{L^2(\Omega)}/||u||_{L^2(\Omega)}');



