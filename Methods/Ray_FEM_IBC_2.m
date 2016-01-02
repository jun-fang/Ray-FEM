function [u,A,v,b,rel_L2_err] = Ray_FEM_IBC_2(node,elem,omega,pde,ray,fquadorder,plt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ray-based FEM for solving the Helmholtz equation with Impedance Boundary Condition(IBC): 
%         -\Delta u - (omega/c)^2 u = f               in D
%         \partial u / \partial n + i omega/c u = g   on \partial D 
%  Version 2: The numbers of rays crossing at each grid node could be different!!! 
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
%   ray: N x 1 cell that contains the ray information
%        stored as the complex form exp(i*ray_angle)
%
%   fquadorder: The order of numerical quadrature
%
%   plt: 1 plot the solution; 0 not plot
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% OUTPUT:
%
%   u: N x 1 vector, stored as the solution value at grid nodes
%   A: Ndof x Ndof assembling matrix 
%   v: Ndof x 1 vector, the solution of Av = b, the coefficients of bases
%   b: Ndof x 1 vector, the right hand side
%   rel_L2_err: relative L2 error
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


N = size(node,1);         % number of grid nodes
NT = size(elem,1);        % number of triangle elements
c = pde.speed(node);      % medium speed 
k = omega./c;             % wavenumber


%% Numerical Quadrature
[lambda,weight] = quadpts(fquadorder);
phi = lambda;             % the value of piesewise linear basis functions at numerical quadrature nodes
nQuad = size(lambda,1);


%% Compute geometric quantities and gradient of local basis
[Dphi,area] = gradbasis(node,elem);


%% Ray information
ray_num = zeros(N,1);     % number of rays at each grid node
ray_dof = zeros(N,1);     % ray_dof(n) = sum(ray_num(1:n)) 

temp = 0;
for n = 1:N
    ray_num(n) = size(ray{n},2);
    ray_dof(n) = temp + ray_num(n);
    temp = ray_dof(n);
end
Ndof = temp;              % degree of freedom

nv = 0;
for nt = 1:NT
    for i = 1:3
        ei = elem(nt,i);
        ni = ray_num(ei);
        for j = 1:3
            ej = elem(nt,j);
            nj = ray_num(ej);
            nv = nv + ni*nj;
        end
    end
end
Nvec = nv*nQuad;          % the length of assembling vectors


%% Assemble the linear system Av = b
b = zeros(Ndof,1);        % right hand side
rows = zeros(Nvec,1);
cols = rows;
vals = rows;
inds = 1;

for nt = 1:NT
    for p = 1:nQuad
        % numerical quadrature nodes in the x-y coordinate
        pxy = lambda(p,1)*node(elem(nt,1),:) ...
            + lambda(p,2)*node(elem(nt,2),:) ...
            + lambda(p,3)*node(elem(nt,3),:);
        fp = pde.f(pxy);
        k2 = (omega/pde.speed(pxy))^2;
        
        for i = 1:3
            ei = elem(nt,i);    % the ei'th grid node corresponding to the i'th vertix of element nt
            ni = ray_num(ei);   % number of rays crossing at grid node ei
            for nii = 1: ni
                gradtempi = - ray{ei}(nii);
                gradtempi = [real(gradtempi), imag(gradtempi)];
                fphasei = gradtempi(1)*pxy(1) + gradtempi(2)*pxy(2);
                ki = k(ei);
                phasei = exp(1i*ki*fphasei);
                
                ii = ray_dof(ei) - ni + nii;
                b(ii) = b(ii) + weight(p)*phi(p,i)*phasei*fp*area(nt);
                
                for j = 1:3
                    ej = elem(nt,j);    % the ej'th grid node corresponding to the j'th vertix of element nt
                    nj = ray_num(ej);   % number of rays crossing at grid node ej
                    
                    for njj = 1: nj
                        gradtempj = ray{ej}(njj);
                        gradtempj = [real(gradtempj), imag(gradtempj)];
                        fphasej = gradtempj(1)*pxy(1) + gradtempj(2)*pxy(2);
                        kj = k(ej);
                        phasej = exp(1i*kj*fphasej);
                        exp_phase = phasei*phasej;
                        
                        jj = ray_dof(ej) - nj + njj;
                  
                        tempA1 = Dphi(nt,1,i)*Dphi(nt,1,j) + Dphi(nt,2,i)*Dphi(nt,2,j);                      
                        
                        tempA2 = 1i*ki*phi(p,i)*(gradtempi(1)*Dphi(nt,1,j) + gradtempi(2)*Dphi(nt,2,j))...
                            + 1i*kj*phi(p,j)*(gradtempj(1)*Dphi(nt,1,i) + gradtempj(2)*Dphi(nt,2,i));
                        
                        tempA3 = ki*kj*phi(p,i)*phi(p,j)...
                            *(gradtempi(1)*gradtempj(1) + gradtempi(2)*gradtempj(2));
                        
                        tempM = k2*phi(p,i)*phi(p,j);
                        
                        Aij = weight(p)*( tempA1 + tempA2 - tempA3 - tempM )*exp_phase*area(nt);
                        
                        rows(inds) = ii;
                        cols(inds) = jj;
                        vals(inds) = Aij;
                        inds = inds + 1;
                        
                    end
                end
            end
        end
    end
end

A = sparse(rows,cols,vals,Ndof,Ndof);


%% Boundaries
[~,bdEdge,~] = findboundary(elem);
Ne = size(bdEdge,1);
el = sqrt(sum((node(bdEdge(:,1),:) - node(bdEdge(:,2),:)).^2,2));

[lambdagN,weightgN] = quadpts1(fquadorder);
phigN = lambdagN;                
nQuadgN = size(lambdagN,1);


%% Modify the linear system Av = b with Impedance Boundary Conditions
for ne = 1:Ne
    for pp = 1:nQuadgN
        ppxy = lambdagN(pp,1)*node(bdEdge(ne,1),:) ...
            + lambdagN(pp,2)*node(bdEdge(ne,2),:);
        gp = pde.g_IBC(ppxy);
        kp = omega/pde.speed(ppxy);
        
        for i = 1:2
            ei = bdEdge(ne,i);    
            ni = ray_num(ei);     
            for nii = 1: ni
                gradtempi = - ray{ei}(nii);
                gradtempi = [real(gradtempi), imag(gradtempi)];
                fphasei = gradtempi(1)*ppxy(1) + gradtempi(2)*ppxy(2);
                ki = k(ei);
                phasei = exp(1i*ki*fphasei);
                
                ii = ray_dof(ei) - ni + nii;
                b(ii) = b(ii) + weightgN(pp)*phigN(pp,i)*phasei*gp*el(ne);   
                
                for j = 1:2
                    ej = bdEdge(ne,j);   
                    nj = ray_num(ej);   
                    
                    for njj = 1: nj
                        gradtempj = ray{ej}(njj);
                        gradtempj = [real(gradtempj), imag(gradtempj)];
                        fphasej = gradtempj(1)*ppxy(1) + gradtempj(2)*ppxy(2);
                        kj = k(ej);
                        phasej = exp(1i*kj*fphasej);
                        exp_phase = phasei*phasej;
                        
                        jj = ray_dof(ej) - nj + njj;
                        Aij = 1i*kp*weightgN(pp)*phigN(pp,i)*phigN(pp,j)*exp_phase*el(ne);
                        A = A + sparse(ii,jj,Aij,Ndof,Ndof);   
                        
                    end
                end
            end
        end
    end
end


%% Solve the linear system Av = b
fprintf('Direct solving time for Av = b: \n');
tic;
v = A\b;
toc;


%% Compute solution values at grid nodes
u = zeros(N,1);
for n = 1:N
    grad = ray{n};
    grad = [real(grad), imag(grad)];    
    tempn = node(n,1)*grad(:,1) + node(n,2)*grad(:,2);
    ni = ray_dof(n) - ray_num(n);
    nii = 1:ray_num(n);
    nii = nii + ni;
    u(n) = sum(v(nii).*exp(1i*k(n)*tempn));
end


%% Plot the solution
if plt
    figure(10);
    uex = pde.ex_u(node);
    showresult(node,elem,real(uex));
    title('Exact Solution');
    
    figure(14);
    showresult(node,elem,real(u));
    title('Ray-based FEM solution with IBC');
end


%% Get relative L2 error
err = zeros(NT,1);
rel_err = err;
for nt = 1:NT
    for p = 1:nQuad
        pxy = lambda(p,1)*node(elem(nt,1),:) ...
            + lambda(p,2)*node(elem(nt,2),:) ...
            + lambda(p,3)*node(elem(nt,3),:);
        temp = 0;
        
        for i = 1:3
            ei = elem(nt,i);     
            ni = ray_num(ei);     
            for nii = 1: ni
                ray_nii = ray{ei}(nii);
                gradtempi = [real(ray_nii), imag(ray_nii)];
                fphasei = gradtempi(1)*pxy(1) + gradtempi(2)*pxy(2);
                ki = k(ei);
                exp_phaseii = exp(1i*ki*fphasei);
                ii = ray_dof(ei) - ni + nii;
                temp = temp + v(ii)*phi(p,i)*exp_phaseii;
            end
        end
        
        abs_err = abs(pde.ex_u(pxy) - temp);
        err(nt) = err(nt) + weight(p)*(abs_err)^2;
        rel_err(nt) = rel_err(nt) + weight(p)*(abs(pde.ex_u(pxy)))^2;
    end
end

err = area.*err;
rel_err = area.*rel_err;
err(isnan(err)) = 0;      % singular values, i.e. uexact(p) = infty, are excluded
err = sqrt(sum(err));
rel_err = sqrt(sum(rel_err));
rel_L2_err = err/rel_err;



