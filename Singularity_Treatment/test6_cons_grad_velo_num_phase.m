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

plt = 0;                   % show solution or not
fquadorder = 3;            % numerical quadrature order
Nray = 1;                  % one ray direction
sec_opt = 0;               % NMLA second order correction or not

high_omega = 100*pi;
low_omega = 2*sqrt(high_omega);


low_wpml = 0.18;
high_wpml = 0.065;

% wavelength
high_wl = 2*pi/high_omega;
low_wl = 2*pi/low_omega;

% mesh size
fh = 1/700;  % fine mesh size
ch = 1/100;                                             % coarse mesh size

sd = 1/2; md = 0.6; ld = 0.9;
Rest = 0.4654; ti = 1;

h = fh; h_c = ch;
tstart = tic;


%% Babich expansion

load('Babich_CGV.mat');

CompressRatio = 5;
Bh = 1/round( 1/(Bh0/CompressRatio) );
a = sd;
Bx = -a: Bh : a;  By = -a: Bh : a;
[BX0, BY0] = meshgrid(Bx0, By0);
[BX, BY] = meshgrid(Bx, By);

% refined phase and amplitude
DD1 = interp2(BX0,BY0,D1,BX,BY,'spline'); % refined amplitude
DD2 = interp2(BX0,BY0,D2,BX,BY,'spline'); % refined amplitude

ttao = interp2(BX0,BY0,tao,BX,BY,'spline'); % refined phase
taox = tao2x ./ (2*tao);   taox(71, 71) = 0;
taoy = tao2y ./ (2*tao);   taoy(71, 71) = 0;
ttaox = interp2(BX0,BY0,taox,BX,BY,'spline'); % refined phase
ttaoy = interp2(BX0,BY0,taoy,BX,BY,'spline'); % refined phase

ray = atan2(ttaoy(:),ttaox(:));
xx = BX(:)-xs;  yy = BY(:)-ys;
rr = sqrt(xx.^2 + yy.^2);
Bray = exp(1i*ray).*(rr>10*eps);
Bray = reshape(Bray, size(BX));


%% Step 1: Solve the Hemholtz equation with the same source but with a relative low frequency sqrt(\omega) by Standard FEM, mesh size \omega*h = constant
fprintf(['-'*ones(1,80) '\n']);
fprintf('Step1: S-FEM, low frequency \n');
tic;
omega = low_omega(ti);              % low frequency
a = ld(ti);                         % large domain
wpml = low_wpml(ti);                % width of PML
sigmaMax = 25/wpml;                 % Maximun absorbtion
[lnode,lelem] = squaremesh([-a,a,-a,a],h);

% smooth part
A = assemble_Helmholtz_matrix_SFEM(lnode,lelem,omega,wpml,sigmaMax,speed,fquadorder);
b = assemble_RHS_SFEM_with_ST(lnode,lelem,xs,ys,omega,wpml,sigmaMax,epsilon,fquadorder);
[~,~,isBdNode] = findboundary(lelem);
freeNode = find(~isBdNode);
lN = size(lnode,1);        u_std = zeros(lN,1);
u_std(freeNode) = A(freeNode,freeNode)\b(freeNode);

% singular part
G1 = 1i*(sqrt(pi)/2)*besselh(0,1,omega*ttao);
G1 = G1.*DD1;

G2 = 1i*(sqrt(pi)/2)*exp(1i*pi)*(2*ttao/omega);
G2 = G2.*besselh(1,1,omega*ttao);
G2 = G2.*DD2;

% x = lnode(:,1); y = lnode(:,2);
% rr = sqrt((x-xs).^2 + (y-ys).^2);
ln = round(sqrt(lN)); 
ub = zeros(ln, ln); 
idx = round((ld-sd)/h)+1:round((ld+sd)/h)+1;
ub(idx, idx) = G1 + G2;
ub = ub(:);
cf = cutoff(epsilon,2*epsilon,lnode,xs,ys);

% low frequency solution: smooth + singularity
u_low = u_std + ub.*cf;
toc;


%% Step 2: Use NMLA to find ray directions d_c with low frequency sqrt(\omega)
fprintf(['-'*ones(1,80) '\n']);
fprintf('Step2: NMLA, low frequency \n');
tic;

% compute numerical derivatives
[ux,uy] = num_derivative(u_low,h,2);

a = md(ti);
[mnode,melem] = squaremesh([-a,a,-a,a],h);
[cnode,celem] = squaremesh([-a,a,-a,a],h_c);
cN = size(cnode,1);
cnumray_angle = zeros(cN,Nray);

% NMLA
for i = 1:cN
    x0 = cnode(i,1);  y0 = cnode(i,2);
    r0 = sqrt((x0-xs)^2 + (y0-ys)^2);
    c0 = speed(cnode(i,:));
    if r0 > (2*epsilon - 3*h_c)
        [cnumray_angle(i)] = NMLA(x0,y0,c0,omega,Rest,lnode,lelem,u_low,ux,uy,[],1/5,Nray,'num',sec_opt,plt);
    else
        cnumray_angle(i) =  ex_ray([x0,y0],xs,ys,0);
    end
end
cnumray = exp(1i*cnumray_angle);
ray = interpolation(cnode,celem,mnode,cnumray);
ray = ray./abs(ray);

% Babich ray directions in the support of the cut-off dunction
x = mnode(:,1); y = mnode(:,2);
rr = sqrt((x-xs).^2 + (y-ys).^2);
mN = length(rr); mn = round(sqrt(mN));
mray = zeros(mn, mn);
idx = round((md-sd)/h)+1:round((md+sd)/h)+1;
mray(idx, idx) = Bray;
mray = mray(:);
ray(rr <= 2*epsilon) = mray(rr <= 2*epsilon);
%     figure(1); ray_field(ray,mnode,10,1/10);
toc;


%% Step 3: Solve the original Helmholtz equation by Ray-based FEM with ray directions d_c
fprintf(['-'*ones(1,80) '\n']);
fprintf('Step3: Ray-FEM, high frequency \n');
tic;
omega = high_omega(ti);
wpml = high_wpml(ti);                % width of PML
sigmaMax = 25/wpml;                 % Maximun absorbtion

% smooth part
option ={ 'Babich_CGV', 'numerical_phase'};
A = assemble_Helmholtz_matrix_RayFEM(mnode,melem,omega,wpml,sigmaMax,speed,ray,fquadorder);
b = assemble_RHS_RayFEM_with_ST(mnode,melem,xs,ys,omega,epsilon,wpml,sigmaMax,ray,speed,fquadorder,option);
uh = RayFEM_direct_solver(mnode,melem,A,b,omega,ray,speed);

% singularity part
G1 = 1i*(sqrt(pi)/2)*besselh(0,1,omega*ttao);
G1 = G1.*DD1;

G2 = 1i*(sqrt(pi)/2)*exp(1i*pi)*(2*ttao/omega);
G2 = G2.*besselh(1,1,omega*ttao);
G2 = G2.*DD2;

ub = zeros(mn, mn); 
idx = round((md-sd)/h)+1:round((md+sd)/h)+1;
ub(idx, idx) = G1 + G2;
ub = ub(:);
cf = cutoff(epsilon,2*epsilon,mnode,xs,ys);

% smooth + singularity
uh1 = uh + ub.*cf;
toc;



%% Step 4: NMLA to find original ray directions d_o with wavenumber k
fprintf(['-'*ones(1,80) '\n']);
fprintf('Step4: NMLA, high frequency \n');
tic;

% compute numerical derivatives
[ux,uy] = num_derivative(uh1,h,2);

a = sd;
[node,elem] = squaremesh([-a,a,-a,a],h);
[cnode,celem] = squaremesh([-a,a,-a,a],h_c);
cN = size(cnode,1);
cnumray_angle = zeros(cN,Nray);

% NMLA
for i = 1:cN
    x0 = cnode(i,1);  y0 = cnode(i,2);
    r0 = sqrt((x0-xs)^2 + (y0-ys)^2);
    c0 = speed(cnode(i,:));
    if r0 > (2*epsilon - 3*h_c)
        [cnumray_angle(i)] = NMLA(x0,y0,c0,omega,Rest,mnode,melem,uh1,ux,uy,[],1/5,Nray,'num',sec_opt,plt);
    else
        cnumray_angle(i) =  ex_ray([x0,y0],xs,ys,0);
    end
end
cnumray = exp(1i*cnumray_angle);
toc;


%% refined ray
h_ray = 1/2400;
[node_ray, elem_ray] = squaremesh([-a,a,-a,a],h_ray);
rray = interpolation(cnode,celem,node_ray,cnumray);
rray = rray./abs(rray);

load('Babich_CGV.mat');

Bh = h_ray;
a = sd;
Bx = -a: Bh : a;  By = -a: Bh : a;
[BX0, BY0] = meshgrid(Bx0, By0);
[BX, BY] = meshgrid(Bx, By);

% refined phase and amplitude
taox = tao2x ./ (2*tao);   taox(71, 71) = 0;
taoy = tao2y ./ (2*tao);   taoy(71, 71) = 0;
ttaox = interp2(BX0,BY0,taox,BX,BY,'spline'); % refined phase
ttaoy = interp2(BX0,BY0,taoy,BX,BY,'spline'); % refined phase

ray = atan2(ttaoy(:),ttaox(:));
xx = BX(:)-xs;  yy = BY(:)-ys;
rr = sqrt(xx.^2 + yy.^2);
Bray = exp(1i*ray).*(rr>10*eps);
% Bray = reshape(Bray, size(BX));
rray(rr <= 2*epsilon) = Bray(rr <= 2*epsilon);
rray = reshape(rray, size(BX));


%% refined Babich
Bh = 1/2800;
a = sd;
Bx = -a: Bh : a;  By = -a: Bh : a;
[BX0, BY0] = meshgrid(Bx0, By0);
[BX, BY] = meshgrid(Bx, By);

% refined phase and amplitude
ttao = interp2(BX0,BY0,tao,BX,BY,'spline'); % refined phase
DD1 = interp2(BX0,BY0,D1,BX,BY,'spline'); % refined amplitude
DD2 = interp2(BX0,BY0,D2,BX,BY,'spline'); % refined amplitude


%% Reference solution
load('C:\Users\Jun Fang\Dropbox\test6_reference_solution.mat');
rh = 1/2800;
[rnode,~] = squaremesh([-a,a,-a,a],rh);
x = rnode(:,1); y = rnode(:,2);
rr = sqrt((x-xs).^2 + (y-ys).^2);


%% h refinement
fh = 1./[400 600 800 1200];

test_num = 4;

% error
max_err = 0*zeros(1,test_num);      % L_inf error of the numerical solution
rel_max_err = 0*zeros(1,test_num);  % relative L_inf error of the numerical solution
l2_err = 0*zeros(1,test_num);       % L_2 error of the numerical solution
rel_l2_err = 0*zeros(1,test_num);   % relative L_2 error of the numerical solution
ref_l2 = 0*zeros(1,test_num);       % reference l2 norm


for ti = 1: test_num
    h = fh(ti);
    [node,elem] = squaremesh([-a,a,-a,a],h);
    N = size(node,1);  n = round(sqrt(N));
    idx = 1:n;
    nh = round(h/h_ray);
    idx = 1 + (idx-1)*nh;
    ray = rray(idx, idx);
    ray = ray(:);
    
   
    %% Step 5: Solve the original Helmholtz equation by Ray-based FEM with ray directions d_o
    fprintf([ '-'*ones(1,80) '\n']);
    fprintf('Step5: Ray-FEM, high frequency \n');
    tic;
    
    % Assembling
    omega = high_omega;
    wpml = 0.1;                         % width of PML
    sigmaMax = 25/wpml;                 % absorption of PML
    
    option ={ 'Babich_CGV', 'numerical_phase'};
    A = assemble_Helmholtz_matrix_RayFEM(node,elem,omega,wpml,sigmaMax,speed,ray,fquadorder);
    b = assemble_RHS_RayFEM_with_ST(node,elem,xs,ys,omega,epsilon,wpml,sigmaMax,ray,speed,fquadorder,option);
    [~,v] = RayFEM_direct_solver(node,elem,A,b,omega,ray,speed);
    toc;
    
    
    %% Compute errors
    fprintf([ '-'*ones(1,80) '\n']);
    fprintf('Compute errors\n');
    tic;
    
    % smooth part
    uh = RayFEM_solution(node,elem,omega,speed,v,ray,rnode);
    
    % singularity part
    G1 = 1i*(sqrt(pi)/2)*besselh(0,1,omega*ttao);
    G1 = G1.*DD1;
    
    G2 = 1i*(sqrt(pi)/2)*exp(1i*pi)*(2*ttao/omega);
    G2 = G2.*besselh(1,1,omega*ttao);
    G2 = G2.*DD2;
    
    ub = G1 + G2;
    ub = ub(:);
    cf = cutoff(epsilon,2*epsilon,rnode,xs,ys);
    
    % smooth + singularity
    uh = uh + ub.*cf;
    
    % Errors
    ur = u_ref;
    du = uh - ur;
    idx = find( ~( (x<=max(x)-wpml).*(x>= min(x)+wpml)...
        .*(y<= max(y)-wpml).*(y>= min(y)+wpml).*(rr>0.1) ) ); % index on PML
    du(idx) = 0;  ur(idx) = 0;
    
    max_err(ti) = norm(du,inf);
    rel_max_err(ti) = norm(du,inf)/norm(ur,inf);
    l2_err(ti) = norm(du)*rh;
    ref_l2(ti) = norm(ur)*rh;
    rel_l2_err(ti) = l2_err(ti)/ref_l2(ti)
    
    toc;

    
end

totaltime = toc(tstart);
fprintf('\n\nTotal running time: % d minutes \n', totaltime/60);

hs = fh;
save('test6_CGV_num_phase.mat', 'omega', 'hs', 'rel_l2_err');

figure(62);
show_convergence_rate(fh(1:test_num),rel_l2_err(1:test_num),'h','Rel L2 err');



