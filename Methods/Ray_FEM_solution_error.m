function [err, rel_L2_err] = Ray_FEM_solution_error(node,elem,omega,pde,ray,v,fquadorder)

%% FEM set up
N = size(node,1);       % number of grid points
NT = size(elem,1);      % number of triangle elements
Nray = size(ray,2);     % number of rays crossing at each grid node
Ndof = N*Nray;          % degree of freedom

c = pde.speed(node);    % medium speed
k = omega./c;           % wavenumber
kk = repmat(k,1,Nray);


%% Numerical Quadrature
[lambda,weight] = quadpts(fquadorder);
phi = lambda;           % the value of piesewise linear basis functions at numerical quadrature nodes
nQuad = size(lambda,1);


%% Compute geometric quantities and gradient of local basis
[Dphi,area] = gradbasis(node,elem);
reparea = repmat(area,Nray,1);


%% Get relative L2 error
err = zeros(NT,1);
rel_err = err;
exp_phaseii = zeros(NT,3,Nray);
for p = 1:nQuad
    pxy = lambda(p,1)*node(elem(:,1),:) ...
        + lambda(p,2)*node(elem(:,2),:) ...
        + lambda(p,3)*node(elem(:,3),:);
    uhp = zeros(NT,1);
    for ni = 1: Nray
        for i = 1:3
            ray_ni = ray(elem(:,i),ni);
            gradtempi = [real(ray_ni), imag(ray_ni)];
            fphasei = gradtempi(:,1).*pxy(:,1) + gradtempi(:,2).*pxy(:,2);
            ki = kk(elem(:,i),ni);
            exp_phaseii(:,i,ni) = exp(1i*ki.*fphasei);
        end
        nii = (ni - 1)*N;
        uhp = uhp + v(elem(:,1) + nii)*phi(p,1).*exp_phaseii(:,1,ni) + ...
            v(elem(:,2) + nii)*phi(p,2).*exp_phaseii(:,2,ni) + ...
            v(elem(:,3) + nii)*phi(p,3).*exp_phaseii(:,3,ni);
    end
    
    abs_err = abs(pde.ex_u(pxy) - uhp);
    err = err + weight(p)*(abs_err).^2;
    rel_err = rel_err + weight(p)*(abs(pde.ex_u(pxy))).^2;
end

err = area.*err;
rel_err = area.*rel_err;
err(isnan(err)) = 0;
err = sqrt(sum(err));
rel_err = sqrt(sum(rel_err));
rel_L2_err = err./rel_err;


