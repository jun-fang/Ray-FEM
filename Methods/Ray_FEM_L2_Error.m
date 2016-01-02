function [rel_L2_err] = Ray_FEM_L2_Error(v,node,elem,omega,c,u_ex,ray,fquadorder)
%% Get the relative L2 error of the Ray-based FEM solution
% u_ex could be a function or reference solution arrary
if isnumeric(u_ex)
    xmax = max(node(:,1));
    xmin = min(node(:,1));
    ymax = max(node(:,2));
    ymin = min(node(:,2));
    refN = length(u_ex);
    refn = round(sqrt(refN));
    refh = (xmax-xmin)/(refn-1);
    [refnode,refelem] = squaremesh([xmin,xmax,ymin,ymax],refh);
end
    
N = size(node,1);
NT = size(elem,1);
Nray = size(ray,2);

% c = pde.speed(node);    % medium speed
k = omega./c;           % wavenumber
kk = repmat(k,1,Nray);

[lambda,weight] = quadpts(fquadorder);
phi = lambda;           
nQuad = size(lambda,1);
[~,area] = gradbasis(node,elem);

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
    if ~isnumeric(u_ex)
        upxy = u_ex(pxy);
    else
        upxy = interpolation(refnode,refelem,pxy,u_ex);
    end       
    abs_err = abs(upxy - uhp);
    err = err + weight(p)*(abs_err).^2;
    rel_err = rel_err + weight(p)*(abs(upxy)).^2;
end

err = area.*err;
rel_err = area.*rel_err;
err(isnan(err)) = 0;
err = sqrt(sum(err));
rel_err = sqrt(sum(rel_err));
rel_L2_err = err./rel_err;
