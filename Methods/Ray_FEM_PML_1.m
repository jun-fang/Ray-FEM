function [u,A,v,b] = Ray_FEM_PML_1(node,elem,omega,wpml,sigmaMax,pde,ray,fquadorder,plt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ray-based FEM for solving the Helmholtz equation with PML: 
%         -\Delta u - (omega/c)^2 u = f               in D
%                                 u = 0               on \partial D 
%  Version 1: Assume the numbers of rays crossing at each grid node are the same!!! 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT:
%
%   node  N x 2 matrix that contains the physical position of each node
%         node(:,1) provides the x coordinate
%         node(:,2) provides the y coordinate
% 
%   elem: NT x 3 matrix that contains the indices of the nodes for each
%         triangle element
%
%   omega: angular frequency
% 
%   wpml: Width of the PML
% 
%   sigmaMax: Maximun absorbtion
%
%   pde: Structure of functions containing the pde information  
%        like speed c(x), right hand side f(x)
%
%   ray: N x Nray matrix that contains the ray information
%        stored as the complex form exp(i*ray_angle)
%
%   fquadorder: The order of numerical quadrature
%
%   plt: 1 plot the solution; 0 not plot
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT:
%
%   u: N x 1 vector, stored as the solution value at grid nodes
%   A: Ndof x Ndof assembling matrix 
%   v: Ndof x 1 vector, the solution of Av = b, the coefficients of bases
%   b: Ndof x 1 vector, the right hand side
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% FEM set up
N = size(node,1);       % number of grid points
NT = size(elem,1);      % number of triangle elements
Nray = size(ray,2);     % number of rays crossing at each grid node
Ndof = N*Nray;          % degree of freedom 

c = pde.speed(node);    % medium speed
k = omega./c;           % wavenumber
kk = repmat(k,1,Nray);

xmax = max(node(:,1));    
xmin = min(node(:,1));
ymax = max(node(:,2));
ymin = min(node(:,2));


%% PML set up
sigmaPML_x = @(x)sigmaMax*( (x-xmin-wpml).^2.*(x < xmin + wpml) + ...
                (x-(xmax-wpml)).^2.*(x > xmax - wpml))/wpml^2;       
sigmaPML_y = @(y) sigmaMax*( (y-ymin-wpml).^2.*(y < ymin + wpml) ...
                + (y-(ymax-wpml)).^2.*(y > ymax - wpml))/wpml^2;
s_x = @(x,y) (1+1i*sigmaPML_y(y)/omega)./(1+1i*sigmaPML_x(x)/omega);       %% s1/s2
s_y = @(x,y) (1+1i*sigmaPML_x(x)/omega)./(1+1i*sigmaPML_y(y)/omega);       %% s2/s1

s_xy = @(x,y) ((1+1i*sigmaPML_x(x)/omega).*(1+1i*sigmaPML_y(y)/omega));    %% 1/(s1*s2)


%% Numerical Quadrature
[lambda,weight] = quadpts(fquadorder);
phi = lambda;           % the value of piesewise linear basis functions at numerical quadrature nodes
nQuad = size(lambda,1);


%% Compute geometric quantities and gradient of local basis
[Dphi,area] = gradbasis(node,elem);
reparea = repmat(area,Nray,1);


%% Assemble the linear system Av = b
Nvec = nQuad*3*3*Nray*Nray*NT;
rows = zeros(Nvec,1);
cols = rows;
vals = rows;
inds = 1:NT;

bt = zeros(NT*Nray,3);
b = zeros(Ndof,1);

fprintf('Assembling time:\n');
tic;
for p = 1:nQuad
    pxy = lambda(p,1)*node(elem(:,1),:) ...
        + lambda(p,2)*node(elem(:,2),:) ...
        + lambda(p,3)*node(elem(:,3),:);
    reppxy = repmat(pxy,Nray,1);
    sx = s_x(pxy(:,1),pxy(:,2));
    sy = s_y(pxy(:,1),pxy(:,2));
    sxy = s_xy(pxy(:,1),pxy(:,2));
    ssxy = repmat(sxy,Nray,1);
    fp = pde.f(reppxy).*ssxy;
    k2 = (omega./pde.speed(pxy)).^2;
    
    for i = 1:3
        gradtempi = - ray(elem(:,i),:);
        gradtempi = gradtempi(:);
        gradtempi = [real(gradtempi), imag(gradtempi)];
        fphasei = gradtempi(:,1).*reppxy(:,1) + gradtempi(:,2).*reppxy(:,2);
        kki = kk(elem(:,i),:);
        kki = kki(:);
        phasei = exp(1i*kki.*fphasei);
        bt(:,i) = bt(:,i) + weight(p)*phi(p,i)*phasei.*fp;
        
        for j = 1:3
            for nii = 1: Nray
                gradtempi = - ray(elem(:,i),nii);
                gradtempi = gradtempi(:);
                gradtempi = [real(gradtempi), imag(gradtempi)];
                fphasei = gradtempi(:,1).*pxy(:,1) + gradtempi(:,2).*pxy(:,2);
                ki = kk(elem(:,i),nii);
                phasei = exp(1i*ki.*fphasei);
                
                for njj = 1: Nray
                    gradtempj = ray(elem(:,j),njj);
                    gradtempj = gradtempj(:);
                    gradtempj = [real(gradtempj), imag(gradtempj)];
                    fphasej = gradtempj(:,1).*pxy(:,1) + gradtempj(:,2).*pxy(:,2);
                    kj = kk(elem(:,j),njj);
                    phasej = exp(1i*kj.*fphasej);
                    exp_phase = phasei.*phasej;
                    
                    tempA1 = sx.*Dphi(:,1,i).*Dphi(:,1,j) + sy.*Dphi(:,2,i).*Dphi(:,2,j); 
                    
                    tempA2 = 1i*ki*phi(p,i).*(sx.*gradtempi(:,1).*Dphi(:,1,j) + sy.*gradtempi(:,2).*Dphi(:,2,j))...
                        + 1i*kj*phi(p,j).*(sx.*gradtempj(:,1).*Dphi(:,1,i) + sy.*gradtempj(:,2).*Dphi(:,2,i));  
                    
                    tempA3 = phi(p,i)*phi(p,j)*ki.*kj.*(sx.*gradtempi(:,1).*gradtempj(:,1)...
                        + sy.*gradtempi(:,2).*gradtempj(:,2));
                    
                    tempM = k2.*sxy*phi(p,i)*phi(p,j);
                    
                    rows(inds) = (nii-1)*N + elem(:,i);
                    cols(inds) = (njj-1)*N + elem(:,j);
                    vals(inds) = weight(p)*(tempA1 + tempA2 - tempA3 - tempM).*exp_phase.*area;
                    inds = inds + NT;
                end
            end
        end
    end
end

A = sparse(rows,cols,vals,Ndof,Ndof);

bt = bt.*repmat(reparea,1,3);
for ni = 1:Nray
    ii = (ni-1)*NT+1:1:ni*NT;
    jj = (ni-1)*N+1:1:ni*N;
    btii = bt(ii,:);
    b(jj) = accumarray(elem(:),btii(:),[N 1]);
end
toc;


%% Boundaries
[~,~,isBdNode] = findboundary(elem);
rep_isBdNode = repmat(isBdNode,1,Nray);
isBdNode = rep_isBdNode(:);
freeNode = find(~isBdNode);


%% Solve the linear system Au = b
v = zeros(Ndof,1);
fprintf('Direct solver time:\n');
tic;
v(freeNode) = A(freeNode,freeNode)\b(freeNode);
toc;


%% Compute solution values at grid nodes
grad = ray(:);
grad = [real(grad),imag(grad)];
repnode = repmat(node,Nray,1);
temp = grad(:,1).*repnode(:,1) + grad(:,2).*repnode(:,2);

u = v.*exp(1i*kk(:).*temp);
u = reshape(u,N,Nray);
u = sum(u,2);


%% Plot the solution
if plt    
    figure(15);
    FJ_showresult(node,elem,real(u));
    title('Ray-based FEM solution with PML');
end


% %% Get relative L2 error
% err = zeros(NT,1);
% rel_err = err;
% exp_phaseii = zeros(NT,3,Nray);
% for p = 1:nQuad
%     pxy = lambda(p,1)*node(elem(:,1),:) ...
%         + lambda(p,2)*node(elem(:,2),:) ...
%         + lambda(p,3)*node(elem(:,3),:);
%     uhp = zeros(NT,1);
%     for ni = 1: Nray
%         for i = 1:3
%             ray_ni = ray(elem(:,i),ni);
%             gradtempi = [real(ray_ni), imag(ray_ni)];
%             fphasei = gradtempi(:,1).*pxy(:,1) + gradtempi(:,2).*pxy(:,2);
%             ki = kk(elem(:,i),ni);
%             exp_phaseii(:,i,ni) = exp(1i*ki.*fphasei);
%         end
%         nii = (ni - 1)*N;
%         uhp = uhp + v(elem(:,1) + nii)*phi(p,1).*exp_phaseii(:,1,ni) + ...
%             v(elem(:,2) + nii)*phi(p,2).*exp_phaseii(:,2,ni) + ...
%             v(elem(:,3) + nii)*phi(p,3).*exp_phaseii(:,3,ni);
%     end
%     
%     abs_err = abs(pde.ex_u(pxy) - uhp);
%     err = err + weight(p)*(abs_err).^2;
%     rel_err = rel_err + weight(p)*(abs(pde.ex_u(pxy))).^2;
% end
% 
% err = area.*err;
% rel_err = area.*rel_err;
% err(isnan(err)) = 0;
% err = sqrt(sum(err));
% rel_err = sqrt(sum(rel_err));
% rel_L2_err = err./rel_err;


