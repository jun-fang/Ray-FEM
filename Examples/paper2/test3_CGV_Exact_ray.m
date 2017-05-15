function [] = test3_CGV_Exact_ray(NPW, test_num)
%% Convergence test for constant gradient of velocity with exact ray

%% Set up
xs = 0; ys = 0;                     % source location
epsilon = 1/(2*pi);                 % cut-off parameter

v0 = 1; G0 = [0.1, -0.2];
speed = @(p) ( v0 + (G0(1)*p(:,1) + G0(2)*p(:,2)) );    % wave speed
speed_min = 0.85;

plt = 0;                   % show solution or not
fquadorder = 3;            % numerical quadrature order
Nray = 1;                  % one ray direction
sec_opt = 0;               % NMLA second order correction or not

high_omega = [100 150 200 300 500 750 1000]*pi;
low_omega = 2*sqrt(high_omega);

% NPW = 4;
% test_num = 3;

% error
max_err = 0*zeros(1,test_num);      % L_inf error of the numerical solution
rel_max_err = 0*zeros(1,test_num);  % relative L_inf error of the numerical solution
l2_err = 0*zeros(1,test_num);       % L_2 error of the numerical solution
rel_l2_err = 0*zeros(1,test_num);   % relative L_2 error of the numerical solution
ref_l2 = 0*zeros(1,test_num);       % reference l2 norm

% mesh size
fh = 1./(10*round(NPW*high_omega/(2*pi*speed_min)/10));      % fine mesh size
ch = 1./(10*round(sqrt(4./fh)/10));                    % coarse mesh size

% width of PML
high_wpml = 0.065*ones(size(high_omega));
low_wpml = 0.185*ones(size(high_omega));


%% Generate the domain sizes
sd = 1/2;
Rest = 0.4654;            % estimate of the distance to the source point

high_r = NMLA_radius(high_omega(1),Rest);
md = sd + high_r + high_wpml;
md = ceil(md*10)/10;      % middle domain size

% Rest = sqrt(2)*md;
low_r = NMLA_radius(low_omega(1),Rest);
ld = md + low_r + low_wpml;
ld = ceil(ld*10)/10;      % large domain size


%% Tests
tstart = tic;
for ti = 1: test_num
    omega = high_omega(ti);
    h = fh(ti);  h_c = ch(ti);
    fprintf(['-'*ones(1,80) '\n']);
    fprintf('\ncase %d: \nomega/(2*pi) = %d,   1/h = %d   1/h_c = %d,  NPW = %d \n',...
        ti, round(omega/(2*pi)), 1/h, 1/h_c, NPW);
    
    
    %% load Babich data
    [~,Bx0,By0,D1,D2,tao,tao2x,tao2y] = load_Babich_data(omega, 'CGV');
    
    a = sd;  Bx = -a: h : a;  By = -a: h : a;
    [BX0, BY0] = meshgrid(Bx0, By0);
    [BX, BY] = meshgrid(Bx, By);
    
    % refined phase and amplitude
    DD1 = interp2(BX0,BY0,D1,BX,BY,'spline'); % refined amplitude
    DD2 = interp2(BX0,BY0,D2,BX,BY,'spline'); % refined amplitude
    
    ttao = interp2(BX0,BY0,tao,BX,BY,'spline'); % refined phase
    mid = round(size(tao,1)/2);
    taox = tao2x ./ (2*tao);   taox(mid, mid) = 0;
    taoy = tao2y ./ (2*tao);   taoy(mid, mid) = 0;
    ttaox = interp2(BX0,BY0,taox,BX,BY,'spline'); % refined phase
    ttaoy = interp2(BX0,BY0,taoy,BX,BY,'spline'); % refined phase
    
    ray = atan2(ttaoy(:),ttaox(:));
    xx = BX(:)-xs;  yy = BY(:)-ys;
    rr = sqrt(xx.^2 + yy.^2);
    Bray = exp(1i*ray).*(rr>10*eps);
    Bray = reshape(Bray, size(BX));
    
    % singularity part
    G1 = 1i*(sqrt(pi)/2)*besselh(0,1,omega*ttao);
    G1 = G1.*DD1;
    
    G2 = 1i*(sqrt(pi)/2)*exp(1i*pi)*(2*ttao/omega);
    G2 = G2.*besselh(1,1,omega*ttao);
    G2 = G2.*DD2;
    
    %% exact ray information
    [node,elem] = squaremesh([-a,a,-a,a],h);
    node0 = [xs, ys];
    x = node(:,1); y = node(:,2);
    rr = sqrt((x-xs).^2 + (y-ys).^2);
    ray = eikonal_cgv(v0, G0, node0, node);
    
    
    %% Step 5: Solve the original Helmholtz equation by Ray-based FEM with ray directions d_o
    fprintf([ '-'*ones(1,80) '\n']);
    fprintf('Step5: Ray-FEM, high frequency \n');
    tic;
    
    omega = high_omega(ti);
    wpml = high_wpml(ti);                % width of PML
    sigmaMax = 25/wpml;                 % Maximun absorbtion
    
    % Assembling
    option ={ 'Babich', 'CGV', 'numerical_phase'};
    A = assemble_Helmholtz_matrix_RayFEM(node,elem,omega,wpml,sigmaMax,speed,ray,fquadorder);
    b = assemble_RHS_RayFEM_with_ST(node,elem,xs,ys,omega,epsilon,wpml,sigmaMax,ray,speed,fquadorder,option);
    uh = RayFEM_direct_solver(node,elem,A,b,omega,ray,speed);
    toc;
    
    
    %% Compute errors
    fprintf([ '-'*ones(1,80) '\n']);
    fprintf('Compute errors\n');
    tic;
    
    % reference solution
    ub = G1 + G2; ub = ub(:);
    cf = cutoff(epsilon,2*epsilon,node,xs,ys);
    ur = (1-cf).*ub;
    ur(rr <= epsilon) = 0;
    
    % Errors
    du = uh - ur;
    idx = find( ~( (x<=max(x)-wpml).*(x>= min(x)+wpml)...
        .*(y<= max(y)-wpml).*(y>= min(y)+wpml) ) ); % index on PML
    du(idx) = 0;  ur(idx) = 0;
    
    max_err(ti) = norm(du,inf);
    rel_max_err(ti) = norm(du,inf)/norm(ur,inf);
    l2_err(ti) = norm(du)*h;
    ref_l2(ti) = norm(ur)*h;
    rel_l2_err(ti) = l2_err(ti)/ref_l2(ti)
    
    toc;
    
    
end

totaltime = toc(tstart);
fprintf('\n\nTotal running time: % d minutes \n', totaltime/60);

nameFile = strcat('results_3_CGV_ExRay_NPW_', num2str(NPW), '.mat');
save(nameFile, 'rel_l2_err', 'NPW', 'high_omega', 'test_num');

% figure(62);
% show_convergence_rate(high_omega(1:test_num),rel_l2_err(1:test_num),'omega','Rel L2 err');



