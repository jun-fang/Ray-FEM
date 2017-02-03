%% Test the convergence rate of Standard FEM

% add path
clear;
addpath(genpath('../../ifem/'));
addpath('../Methods/');
addpath('../Helmholtz_data/');
addpath('../Plots_Prints/');

% set up
pde = Helmholtz_data_one_point_source;    % pde data
xs = 2; ys = 2;       % source location
omega = 10*pi;     % frequency
h = 1/50;                % mesh size
a = 1/2;                  % domain size [-a,a]^2
fquadorder = 3;    % the order of numerical quadrature
solver = 'DIR';      % the type of solver for the resulting linear system
plt = 0;                   % plot the solution or not


% tests
nt = 4;                         % number of tests
hs = zeros(nt,1);          % record mesh size
errors = zeros(nt,1);    % record error

for ti = 1:nt
    h = h/2; hs(ti) = h;
    [node, elem] = squaremesh([-a,a,-a,a], h);
    
    [~,~,~,rel_L2_err] = Standard_FEM_IBC(node,elem,omega,pde,fquadorder,solver,plt,xs,ys);
    errors(ti) = rel_L2_err;
end

% plot
figure(10);
show_convergence_rate(hs, errors);