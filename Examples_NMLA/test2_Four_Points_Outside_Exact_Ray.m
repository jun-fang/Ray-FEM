%% Test One Point Source Problem (outside domain): Iterative Idea
% uexact = sqrt(k)*besselh(0,1,k*sqrt((x-2)^2 + (y-2)^2))
clear;

%% Load source data
global xs ys omega a;
xs = 20; ys = 20;
pde = Helmholtz_data_2;

%% Set up
plt = 0;                   % show solution or not
fquadorder = 9;            % numerical quadrature order
solver = 'DIR';            % linear system solver
pct = 1/8;
Nray = 4;
sec_opt = 0;               % NMLA second order correction or not

rec_N = 8;                 % we test rec_N examples

% record h and omega
rec_omega = zeros(1,rec_N);
rec_h = rec_omega;

% record the error of S-FEM for low frequency equation
rec_low_err = rec_omega;

% record the L2 error of the numerical angle estimation
rec_ang_err1 = rec_omega;
rec_ang_err2 = rec_omega;

% record the error of Numerical Ray-based FEM
rec_NR_err1 = rec_omega;
rec_NR_err2 = rec_omega;

% record the error of Exact Ray-based FEM
rec_ER_err = rec_omega;


% lg_a,md_a,sm_a have to be chosen to match the mesh size h and ch
% (for example md_a/ch = integer), and the real radius in NMLA as well
% (r1+md_a < lg_a, r2 + sm_a < md_a)
% lg_a = 2.1;  md_a = 1.1;  sm_a = 0.5;
% Rest1 = 20;  Rest2 = 20;
NPW = 6;

% cp_omega = [10 20 40 60 80 120 160]*pi;
% cp_omega = [10 20 40 80]*pi;
cp_omega = [10 20 30 40 60 80 120 160]*pi;

sm_a = 1/2;



tstart = tic;
for rec_i = 1: rec_N
    high_omega = cp_omega(rec_i);
    low_omega = sqrt(high_omega);
    h = 1/(NPW*round(high_omega/(2*pi)));
    ch = 1/(10*round(low_omega/(2*pi)));
    

    fprintf(['-'*ones(1,80) '\n']);
    fprintf(['-'*ones(1,80) '\n']);
    fprintf('\ncase %d: \nomega/(2*pi) = %d,   1/h = %d   1/ch = %d,  NPW = %d \n\n',...
        rec_i, round(high_omega/(2*pi)), round(1/h), round(1/ch), NPW);
   
    rec_omega(rec_i) = high_omega;
    rec_h(rec_i) = h;
    
    
    %% Exact Ray-based FEM:
    if (1)
        fprintf('\nExact Ray-based FEM: \n');
        tic;
        omega = high_omega;
        a = sm_a;
        [node,elem] = squaremesh([-a,a,-a,a],h);
        ray = pde.ray(node);
        uh = Ray_FEM_IBC_1(node,elem,omega,pde,ray,fquadorder,plt);
        exu = pde.ex_u(node);
        rec_ER_err(rec_i) = h*norm(uh-exu,2);
        toc;
    end
    
    
end

totaltime = toc(tstart);
fprintf('\n\nTotal running time: % d minutes \n', totaltime/60);



%% record and print results
rec = [rec_omega;rec_h;rec_low_err;rec_ER_err;rec_ang_err1;rec_ang_err2;rec_NR_err1;rec_NR_err2];

fprintf( ['\n\n' '-'*ones(1,80) '\n']);
fprintf( 'Exact Ray-based FEM:');
fprintf( '\n\nRelative L2 error:      ');
fprintf( '&  %1.2d  ',rec_ER_err);
fprintf( ['\n' '-'*ones(1,80) '\n']);

save('result2_ex2.mat','omega','rayerr1','rayerr2','uerr1','uerr2','uerr_exray');

