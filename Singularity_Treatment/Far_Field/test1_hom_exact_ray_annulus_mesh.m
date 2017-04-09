%% Compute errors of Ray-FEM solution in homogeneous medium

clear;

xs = 0;   ys = 0;                  % point source location
speed = @(x) ones(size(x,1),1);    % medium speed

fquadorder = 3;                    % numerical quadrature order

epsilon = 50/(80*pi);    

nt = 4;
omegas = [120 160 240 320]*pi;

rel_l2_err = zeros(1,nt);           % record relative l2 errors


for ni = 1:nt
    
    omega = omegas(ni);
    
    wavelength = 2*pi/omega;
    NPW = 4;
    h = wavelength/NPW;
    h = 1/round(1/h);
    
    wpml = 0.1;%3*round(wavelength/h)*h;                 % width of PML
    sigmaMax = 25/wpml;                % Maximun absorbtion
    
    a = 0.5; bd = 0.1;
    % [node,elem] = squaremesh([-a,a,-a,a],h);
    [node,elem] = squaremesh_annulus([-a,a,-a,a],[-bd,bd,-bd,bd],h);
    Nray = 1; N = size(node,1); Ndof = N;
    
    %% Exact ray information
    x = node(:,1);  y = node(:,2);
    rr = sqrt((x-xs).^2 + (y-ys).^2);
    ray = atan2(y-ys, x-xs);
    ray = exp(1i*ray).*(rr>10*eps);
    
    tic;
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
figure(22);
show_convergence_rate(omegas(1:nt), rel_l2_err,'omega','||u - u_h||_{L^2(\Omega)}/||u||_{L^2(\Omega)}');




