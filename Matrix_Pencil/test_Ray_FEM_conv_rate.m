global omega xs ys a;
xs = 100; ys = 100;                 % source location
aa = 1/2;                         % computational domain [-1/2,1/2]^2

pde = MPM_data5;                  % one point source pde structure
Nray = 2;                         % one ray direction at every grid point
plt = 0;                          % plot solution or not
fquadorder = 6;                   % numerical quadrature order
solver = 'DIR';                   % direct method to solve Au=b by u=A\b  

%% record frequencies and errors
rec_omega = [10 20 30 40 50 60 70 80]*pi;   
rec_err9 = zeros(length(rec_omega),1);
rec_err10 = zeros(length(rec_omega),1);
rec_err11 = zeros(length(rec_omega),1);

nn = 4;

for j = 1:nn%length(rec_omega)
    tic;
    j
    low_omega = sqrt(rec_omega(j));    % low frequency
    omega = rec_omega(j);              % high frequency
    NPW = 8;                           % number of grid points per wavelength  
    h = 1/round((NPW*omega)/(2*pi));   % mesh size
    a = 1/2;
    [node,elem] = squaremesh([-a,a,-a,a],h);
    exray = pde.ray(node);
    
    ray = exp(1i*exray);
    [uh,~,~,~,rel_L2_err] = Ray_FEM_IBC_1(node,elem,omega,pde,ray,fquadorder,plt);
    ue = pde.ex_u(node);
    rec_err9(j) = norm(ue-uh,inf);        % L^{infty} error for high frequency solution
    rec_err10(j) = rel_L2_err; 
    rec_err11(j) = norm(uh-ue)/norm(ue);
    toc;
end

subplot(1,3,1);
showrate(rec_omega(1:nn),rec_err9(1:nn)');
xlabel('\omega');
ylabel('ER-FEM Err_{u, {L^{2}(\Omega)}}');


subplot(1,3,2);
showrate(rec_omega(1:nn),rec_err10(1:nn)');
xlabel('\omega');
ylabel('ER-FEM Err_{u, {L^{\infty}(\Omega)}}');

subplot(1,3,3);
showrate(rec_omega(1:nn),rec_err11(1:nn)');
xlabel('\omega');
ylabel('ER-FEM Err_{u, {L^{\infty}(\Omega)}}');