%% Convergence test for gravity Helmholtz equation with exact ray information


addpath(genpath('../../ifem/'));
addpath('../Methods/');
addpath('../NMLA/');
addpath('../Cutoff_Functions/')
addpath('../Plots_Prints/');


xs = 0; ys = 0;                     % source location 
epsilon = 50/(80*pi);               % cut-off parameter   
omega0 = 80*pi;
E = omega0;
speed = @(p) omega0./sqrt(E + (p(:,1)));    % wave speed

h = 1/16;
[node,elem] = squaremesh([-1,1,-1,1],h);
c = speed(node);



T = eikonal_cgss(S0, grad0, node0, node);

showsolution(node,elem,c);





wpml = 0.1;                         % width of PML  
sigmaMax = 25/wpml;                 % absorption of PML  
fquadorder = 3;                     % numerical quadrature order 
a = 1/2;                            % computational domain [-a,a]^2


nt = 3;                             % number of tests
errors = zeros(1,nt);
rhss = zeros(1,nt);
omegas = pi*[120,160,240,320];      % omega's
NPW = 4;                            % grid number per wavelength



for ii = 1:nt
    ii
    omega = omegas(ii);
    wl = 2*pi/omega;
    h = 1/round(1/(wl/NPW));
    1/h
    [node,elem] = squaremesh([-a,a,-a,a],h);
    
    %% Exact ray information
    xx = node(:,1)-xs;  yy = node(:,2)-ys;
    rr = sqrt(xx.^2 + yy.^2);
    ray = atan2(yy,xx);
    ray = exp(1i*ray).*(rr>10*eps);
    
    %% right hand side
    rhs = sing_rhs_homo(epsilon,omega,node,xs,ys);
    rhss(ii) = norm(rhs)*h;
    
    
    tic;
    % Ray-FEM solution with singularity treatment 
    [u,~,~,v] = RayFEM_singularity(node,elem,xs,ys,omega,epsilon,wpml,sigmaMax,ray,speed,@sing_rhs_homo,fquadorder);
    
    % get the exact solution
    x = node(:,1); y = node(:,2);
    rr = sqrt((x-xs).^2 + (y-ys).^2);
    
    ub = 1i/4*besselh(0,1,omega*rr);
    cf = cutoff(epsilon,2*epsilon,node,xs,ys);
    uex = (1-cf).*ub;
    uex(rr<epsilon) = 0;
    du = u - uex;
    
    idx = find( (x<=max(x)-wpml).*(x>= min(x)+wpml)...
        .*(y<= max(y)-wpml).*(y>= min(y)+wpml) );
    du_phy = du(idx);
    
    dd = 0*du; dd(idx) = du(idx);
    
    % compute the error
    [err, rel_L2_err] = RayFEM_smooth_solution_error(node,elem,xs,ys,omega,epsilon,wpml,ray,speed,v,9);
    
    errors(ii) = err;%norm(du_phy)*h;%/norm(uex(idx));
    toc;
end


%% plot
figure(22);
subplot(1,2,1);
show_convergence_rate(omegas(1:nt),rhss,'omega','||f||_{L^2(\Omega)}');
subplot(1,2,2);
show_convergence_rate(omegas(1:nt),errors,'omega','||u - u_h||_{L^2(\Omega)}');



