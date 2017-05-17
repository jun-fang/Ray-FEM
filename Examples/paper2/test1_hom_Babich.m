%% Convergence test for constant gradient of velocity with exact ray

%% Set up
xs = 0; ys = 0;                     % source location
epsilon = 1/(2*pi);                 % cut-off parameter

v0 = 1; G0 = [0.1, -0.2];
speed = @(p) ones(size(p(:,1)));    % wave speed

plt = 0;                   % show solution or not
fquadorder = 3;            % numerical quadrature order
Nray = 1;                  % one ray direction
sec_opt = 0;               % NMLA second order correction or not
sd = 1/2;

high_omega = [250 400 600 1000 1500]*pi;

test_num = 5;

% error
max_err = 0*zeros(1,test_num);      % L_inf error of the numerical solution
rel_max_err = 0*zeros(1,test_num);  % relative L_inf error of the numerical solution
l2_err = 0*zeros(1,test_num);       % L_2 error of the numerical solution
rel_l2_err = 0*zeros(1,test_num);   % relative L_2 error of the numerical solution
ref_l2 = 0*zeros(1,test_num);       % reference l2 norm


%% Tests
tstart = tic;
for ti = 1: test_num
    
    omega = high_omega(ti);
    
    %% load Babich data
    [Bh0,Bx0,By0,D1,D2,tao,tao2x,tao2y] = load_Babich_data(omega, 'Homo');
    
    [BX0, BY0] = meshgrid(Bx0, By0);
    rr = sqrt((BX0(:)-xs).^2 + (BY0(:)-ys).^2);    
    uex = 1i/4*besselh(0,1,omega*rr);
    
    % singularity part
    G1 = 1i*(sqrt(pi)/2)*besselh(0,1,omega*tao);
    G1 = G1.*D1;
    
    G2 = 1i*(sqrt(pi)/2)*exp(1i*pi)*(2*tao/omega);
    G2 = G2.*besselh(1,1,omega*tao);
    G2 = G2.*D2;
    
    ub = G1 + G2;
    ub = ub(:);
    
    du = ub - uex;
    
    du = tao(:) - rr; uex = rr;
    du(rr < 0.1) = 0; du(rr > 0.3) = 0;
    uex(rr < 0.1) = 0; uex(rr > 0.3) = 0;
    
    max_err(ti) = norm(du,inf);
    rel_max_err(ti) = norm(du,inf)/norm(uex,inf);
    l2_err(ti) = norm(du)*Bh0;
    ref_l2(ti) = norm(uex)*Bh0;
    rel_l2_err(ti) = l2_err(ti)/ref_l2(ti)
        
    
end

totaltime = toc(tstart);
fprintf('\n\nTotal running time: % d minutes \n', totaltime/60);


figure(11);
show_convergence_rate(high_omega(1:test_num),rel_l2_err(1:test_num),'omega','Rel L2 err');

figure(12);
show_convergence_rate(high_omega(1:test_num),max_err(1:test_num),'omega','Rel L2 err');




