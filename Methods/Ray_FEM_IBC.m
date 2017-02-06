function [u,A,v,b,rel_L2_err,err] = Ray_FEM_IBC(node,elem,omega,pde,ray,fquadorder,plt,xs,ys)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ray-based FEM for solving the Helmholtz equation with Impedance Boundary Condition(IBC):
%         -\Delta u - (omega/c)^2 u = f               in D
%         \partial u / \partial n + i omega/c u = g   on \partial D
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
%   xs,ys: source location if it is point source problem
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT:
%
%   u: N x 1 vector, stored as the solution value at grid nodes
%   A: Ndof x Ndof assembling matrix
%   v: Ndof x 1 vector, the solution of Av = b, the coefficients of bases
%   b: Ndof x 1 vector, the right hand side
%   rel_L2_err: relative L2 error
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% FEM set up
N = size(node,1);       % number of grid points
NT = size(elem,1);      % number of triangle elements
Nray = size(ray,2);     % number of rays crossing at each grid node
Ndof = N*Nray;          % degree of freedom
a = max(node(:,1));

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


%% Assemble the linear system Av = b

threshold = 50*pi;
if fquadorder < 6 || omega <= threshold % fast but more memory
    Nvec = nQuad*3*3*Nray*Nray*NT;
    rows = zeros(Nvec,1);
    cols = rows;
    vals = rows;
    inds = 1:NT;
end

A = sparse(Ndof,Ndof);

bt = zeros(NT*Nray,3);
b = zeros(Ndof,1);

% fprintf('Assembling time:\n');
% tic;
for p = 1:nQuad
    pxy = lambda(p,1)*node(elem(:,1),:) ...
        + lambda(p,2)*node(elem(:,2),:) ...
        + lambda(p,3)*node(elem(:,3),:);
    reppxy = repmat(pxy,Nray,1);
    fp = pde.f(reppxy,xs,ys);
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
                    
                    tempA1 = Dphi(:,1,i).*Dphi(:,1,j) + Dphi(:,2,i).*Dphi(:,2,j);
                    
                    tempA2 = 1i*ki*phi(p,i).*(gradtempi(:,1).*Dphi(:,1,j) + gradtempi(:,2).*Dphi(:,2,j))...
                        + 1i*kj*phi(p,j).*(gradtempj(:,1).*Dphi(:,1,i) + gradtempj(:,2).*Dphi(:,2,i));
                    
                    tempA3 = phi(p,i)*phi(p,j)*ki.*kj.*(gradtempi(:,1).*gradtempj(:,1)...
                        + gradtempi(:,2).*gradtempj(:,2));
                    
                    tempM = k2*phi(p,i)*phi(p,j);
                    
                    if fquadorder >= 6 && omega > threshold % slow but less memory
                        ii = (nii-1)*N + elem(:,i);
                        jj = (njj-1)*N + elem(:,j);
                        Aij = weight(p)*(tempA1 + tempA2 - tempA3 - tempM).*exp_phase.*area;
                        A = A + sparse(ii,jj,Aij,Ndof,Ndof);
                    end
                    
                    if fquadorder < 6 || omega <= threshold % fast but more memory
                        rows(inds) = (nii-1)*N + elem(:,i);
                        cols(inds) = (njj-1)*N + elem(:,j);
                        vals(inds) = weight(p)*(tempA1 + tempA2 - tempA3 - tempM).*exp_phase.*area;
                        inds = inds + NT;
                    end
                end
            end
        end
    end
end

if fquadorder < 6 || omega <= threshold % fast but more memory
    A = sparse(rows,cols,vals,Ndof,Ndof);
end

clear Aij Dphi fphasei fphasej phasei phasej exp_phase gradtempi gradtempj;
clear k2 ki kj kki pxy fp reppxy tempA1 tempA2 tempA3 tempM;

% u = 0;v=0;rel_L2_err=0;
% return;

bt = bt.*repmat(reparea,1,3);
for ni = 1:Nray
    ii = (ni-1)*NT+1:1:ni*NT;
    jj = (ni-1)*N+1:1:ni*N;
    btii = bt(ii,:);
    b(jj) = accumarray(elem(:),btii(:),[N 1]);
end

clear bt ii jj bt btii reparea;

% toc;


%% Boundaries
[~,bdEdge,~] = findboundary(elem);
Ne = size(bdEdge,1);      % number of boundary edges
el = sqrt(sum((node(bdEdge(:,1),:) - node(bdEdge(:,2),:)).^2,2));
repel = repmat(el,Nray,1);

[lambdagN,weightgN] = quadpts1(fquadorder);
phigN = lambdagN;
nQuadgN = size(lambdagN,1);


%% Modify the linear system Av = b with Impedance Boundary Condition
rows = zeros(nQuadgN*2*2*Nray*Nray*Ne, 1);
cols = rows;
vals = rows;
inds = 1:Ne;

ge = zeros(Ne*Nray,2);

for pp = 1:nQuadgN
    ppxy = lambdagN(pp,1)*node(bdEdge(:,1),:) ...
        + lambdagN(pp,2)*node(bdEdge(:,2),:);
    kxy = omega./pde.speed(ppxy);
    reppxy = repmat(ppxy,Nray,1);
    gp = pde.g_IBC(reppxy,xs,ys, omega,a);
    
    for i = 1:2
        gradtempi = - ray(bdEdge(:,i),:);
        gradtempi = gradtempi(:);
        gradtempi = [real(gradtempi), imag(gradtempi)];
        fphasei = gradtempi(:,1).*reppxy(:,1) + gradtempi(:,2).*reppxy(:,2);
        kki = kk(bdEdge(:,i),:);
        kki = kki(:);
        phasei = exp(1i*kki.*fphasei);
        ge(:,i) = ge(:,i) + weightgN(pp)*phigN(pp,i)*phasei.*gp;
        
        for j = 1:2
            for nii = 1:Nray
                gradtempi = - ray(bdEdge(:,i),nii);
                gradtempi = gradtempi(:);
                gradtempi = [real(gradtempi), imag(gradtempi)];
                fphasei = gradtempi(:,1).*ppxy(:,1) + gradtempi(:,2).*ppxy(:,2);
                ki = kk(bdEdge(:,i),nii);
                phasei = exp(1i*ki.*fphasei);
                
                for njj = 1: Nray
                    gradtempj = ray(bdEdge(:,j),njj);
                    gradtempj = gradtempj(:);
                    gradtempj = [real(gradtempj), imag(gradtempj)];
                    fphasej = gradtempj(:,1).*ppxy(:,1) + gradtempj(:,2).*ppxy(:,2);
                    kj = kk(bdEdge(:,j),njj);
                    phasej = exp(1i*kj.*fphasej);
                    exp_phase = phasei.*phasej;
                    
                    rows(inds) = (nii-1)*N + bdEdge(:,i);
                    cols(inds) = (njj-1)*N + bdEdge(:,j);
                    vals(inds) = 1i*kxy*weightgN(pp)*phigN(pp,i)*phigN(pp,j).*exp_phase.*el;
                    inds = inds + Ne;
                end
            end
        end
    end
end

A = A + sparse(rows,cols,vals,Ndof,Ndof);

ge = ge.*repmat(repel,1,2);
for ni = 1:Nray
    ii = (ni-1)*Ne+1:1:ni*Ne;
    jj = (ni-1)*N+1:1:ni*N;
    geii = ge(ii,:);
    b(jj) = b(jj) + accumarray(bdEdge(:), geii(:),[N 1]);
end


%% Solve the linear system Av = b
% fprintf('Direct solving time for Av = b: \n');
% tic;
v = A\b;
% toc;


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
    figure(10);
    uex = pde.u_ex(node,xs,ys,omega);
    showresult(node,elem,real(uex));
    title('Exact Solution');
    
    figure(13);
    showresult(node,elem,real(u));
    title('Ray-based FEM with IBC');
end


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
    
    abs_err = abs(pde.u_ex(pxy,xs,ys,omega) - uhp);
    err = err + weight(p)*(abs_err).^2;
    rel_err = rel_err + weight(p)*(abs(pde.u_ex(pxy,xs,ys,omega))).^2;
end

err = area.*err;
rel_err = area.*rel_err;
err(isnan(err)) = 0;
err = sqrt(sum(err));
rel_err = sqrt(sum(rel_err));
rel_L2_err = err./rel_err;


