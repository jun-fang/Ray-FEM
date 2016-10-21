function [u,A,v,b,rel_L2_err] = Phase_FEM_IBC(node,elem,omega,pde,fquadorder,plt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Phase-based FEM for solving the Helmholtz equation with Impedance Boundary Condition(IBC): 
%         -\Delta u - (omega/c)^2 u = f               in D
%         \partial u / \partial n + i omega/c u = g   on \partial D 
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
%        like speed c(x), right hand side f(x), exact phase
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
%   rel_L2_err: relative L2 error
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% FEM set up
N = size(node,1);         % number of grid points
NT = size(elem,1);        % number of triangle elements
Nphase = pde.Nphase;      % number of phase
Ndof = N*Nphase;          % degree of freedom 


%% Numerical quadrature
[lambda,weight] = quadpts(fquadorder);
phi = lambda;                
nQuad = size(lambda,1);


%% Compute geometric quantities and gradient of local basis
[Dphi,area] = gradbasis(node,elem);
reparea = repmat(area,Nphase,1);


%% Assemble the linear system Av = b
Nvec = nQuad*3*3*Nphase*Nphase*NT;
rows = zeros(Nvec,1);
cols = rows;
vals = rows;
inds = 1:NT;

bt = zeros(NT*Nphase,3);
b = zeros(Ndof,1);

fprintf('Assembling time: \n');
tic;
for p = 1:nQuad
    pxy = lambda(p,1)*node(elem(:,1),:) ...
        + lambda(p,2)*node(elem(:,2),:) ...
        + lambda(p,3)*node(elem(:,3),:);
    reppxy = repmat(pxy,Nphase,1);
    k = omega./pde.speed(pxy);
    k2 = k.^2;
    kk = repmat(k,1,Nphase);
    fp = pde.f(reppxy);
    phase = pde.phase(pxy);           % NT x Nphase matrix
    gradphase = pde.gradphase(pxy);   % NT x Nphase matrix, stored as complex form
    
    for i = 1:3
        phasei = - phase;
        exphasei = exp(1i*kk.*phasei);
        exphasei = exphasei(:); 
        bt(:,i) = bt(:,i) + weight(p)*phi(p,i)*exphasei.*fp;
        
        for j = 1:3
            for nii = 1: Nphase
                phasei = - phase(:,nii);
                exphasei = exp(1i*k.*phasei);
                
                gradphasei = - gradphase(:,nii);
                gradphasei = [real(gradphasei), imag(gradphasei)];
                                               
                for njj = 1: Nphase
                    phasej = phase(:,njj);
                    exphasej = exp(1i*k.*phasej);
                    
                    gradphasej = gradphase(:,njj);
                    gradphasej = [real(gradphasej), imag(gradphasej)];
                    
                    exp_phase = exphasei.*exphasej;
                    
                    tempA1 = Dphi(:,1,i).*Dphi(:,1,j) + Dphi(:,2,i).*Dphi(:,2,j); 
                    
                    tempA2 = 1i*k*phi(p,i).*(gradphasei(:,1).*Dphi(:,1,j) + gradphasei(:,2).*Dphi(:,2,j))...
                        + 1i*k*phi(p,j).*(gradphasej(:,1).*Dphi(:,1,i) + gradphasej(:,2).*Dphi(:,2,i));  
                    
                    tempA3 = phi(p,i)*phi(p,j)*k2.*(gradphasei(:,1).*gradphasej(:,1)...
                        + gradphasei(:,2).*gradphasej(:,2));
                    
                    tempM = phi(p,i)*phi(p,j)*k2;
                    
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
for ni = 1:Nphase
    ii = (ni-1)*NT+1:1:ni*NT;
    jj = (ni-1)*N+1:1:ni*N;
    btii = bt(ii,:);
    b(jj) = accumarray(elem(:),btii(:),[N 1]);
end
toc;


%% Boundaries
[~,bdEdge,~] = findboundary(elem);
Ne = size(bdEdge,1);      % number of boundary edges
el = sqrt(sum((node(bdEdge(:,1),:) - node(bdEdge(:,2),:)).^2,2));
repel = repmat(el,Nphase,1);

[lambdagN,weightgN] = quadpts1(fquadorder);
phigN = lambdagN;                
nQuadgN = size(lambdagN,1);


%% Modify the linear system Av = b with Impedance Boundary Condition
rows = zeros(nQuadgN*2*2*Nphase*Nphase*Ne, 1);
cols = rows;
vals = rows;
inds = 1:Ne;
ge = zeros(Ne*Nphase,2);

for pp = 1:nQuadgN
    ppxy = lambdagN(pp,1)*node(bdEdge(:,1),:) ...
        + lambdagN(pp,2)*node(bdEdge(:,2),:);
    reppxy = repmat(ppxy,Nphase,1);
    gp = pde.g_IBC(reppxy);
    k = omega./pde.speed(ppxy);
    kk = repmat(k,1,Nphase);
    phase = pde.phase(ppxy);           % Ne x Nphase matrix
    
    for i = 1:2
        phasei = - phase;
        exphasei = exp(1i*kk.*phasei);     
        exphasei = exphasei(:); 
        ge(:,i) = ge(:,i) + weightgN(pp)*phigN(pp,i)*exphasei.*gp;
        
        for j = 1:2
            for nii = 1:Nphase
                phasei = - phase(:,nii);
                exphasei = exp(1i*k.*phasei);             
                for njj = 1: Nphase
                    phasej = phase(:,njj);
                    exphasej = exp(1i*k.*phasej);
                    exp_phase = exphasei.*exphasej;
                    
                    rows(inds) = (nii-1)*N + bdEdge(:,i);
                    cols(inds) = (njj-1)*N + bdEdge(:,j);
                    vals(inds) = 1i*k*weightgN(pp)*phigN(pp,i)*phigN(pp,j).*exp_phase.*el;
                    inds = inds + Ne;
                end
            end
        end
    end
end

A = A + sparse(rows,cols,vals,Ndof,Ndof);

ge = ge.*repmat(repel,1,2);
for ni = 1:Nphase
    ii = (ni-1)*Ne+1:1:ni*Ne;
    jj = (ni-1)*N+1:1:ni*N;
    geii = ge(ii,:);
    b(jj) = b(jj) + accumarray(bdEdge(:), geii(:),[N 1]);
end


%% Solve the linear system Av = b
fprintf('Direct solving time for Av = b: \n');
tic;
v = A\b;
toc;


%% Compute solution values at grid nodes
phase = pde.phase(node);
k = omega./pde.speed(node);
kk = repmat(k,1,Nphase);
exp_phase = exp(1i*kk.*phase);

u = v.*exp_phase(:);
u = reshape(u,N,Nphase);
u = sum(u,2);


%% Plot the solution
if plt
    figure(10);
    uex = pde.ex_u(node);
    showresult(node,elem,real(uex));
    title('Exact Solution');
    
    figure(16);
    showresult(node,elem,real(u));
    title('Phase-based FEM solution with IBC');
end


%% Get relative L2 error
err = zeros(NT,1);
rel_err = err;
exp_phaseii = zeros(NT,3,Nphase);
for p = 1:nQuad
    pxy = lambda(p,1)*node(elem(:,1),:) ...
        + lambda(p,2)*node(elem(:,2),:) ...
        + lambda(p,3)*node(elem(:,3),:);
    k = omega./pde.speed(pxy);
    phase = pde.phase(pxy);           % Ne x Nphase matrix
    uhp = zeros(NT,1);
    for ni = 1: Nphase
        for i = 1:3
            phasei = phase(:,ni);
            exp_phaseii(:,i,ni) = exp(1i*k.*phasei);
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
