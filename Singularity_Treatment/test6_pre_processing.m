%% tes6 reference solution with analytical phase

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


omega = 100*pi;
wpml = 0.1;                         % width of PML
sigmaMax = 25/wpml;                 % absorption of PML
fquadorder = 3;                     % numerical quadrature order
a = 1/2;                            % computational domain [-a,a]^2



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

%% reference solution
[node, elem] = squaremesh([-a, a, -a, a], Bh);

% Assembling
%tic;
option ={ 'Babich_CGV', 'exact_phase'};
tic;
A = assemble_Helmholtz_matrix_RayFEM(node,elem,omega,wpml,sigmaMax,speed,ray,fquadorder);
toc;
tic;
b = assemble_RHS_RayFEM_with_ST(node,elem,xs,ys,omega,epsilon,wpml,sigmaMax,ray,speed,fquadorder,option);
toc;
tic;
uh = RayFEM_direct_solver(node,elem,A,b,omega,ray,speed);
toc;

clear A b elem;

cf = cutoff(epsilon,2*epsilon,node,xs,ys);
u_ref = uh + ub.*cf;
h_ref = Bh;
%showsolution(node, elem, real(u), 2); colorbar;
save('test6_reference_solution.mat', 'u_ref', 'h_ref');
