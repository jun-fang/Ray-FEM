function [u,A,v,b] = Ray_FEM_PML_2(node,elem,omega,wpml,sigmaMax,pde,ray,fquadorder,plt)
%% Ray-based FEM for solving Helmholtz equation with PML: 
%         -\Delta u - (omega/c)^2 u = f               in D
%                                 u = 0               on \partial D 
%  Version 2: The numbers of rays crossing at each grid node could be different!!! 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT:
%
%   node: N x 2 matrix that contains the physical position of each node
%         node(:,1) provides the x coordinate
%         node(:,2) provides the y coordinate
%
%   elem: NT x 3 matrix that contains the indices of the nodes for each
%         triangle element
%
%   omega: Angular frequency
% 
%   wpml: Width of the PML
% 
%   sigmaMax: Maximun absorbtion
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
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% FEM set up
N = size(node,1);         % number of grid nodes
NT = size(elem,1);        % number of triangle elements
c = pde.speed(node);      % medium speed 
k = omega./c;             % wavenumber

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
phi = lambda;             % the value of piesewise linear basis functions at numerical quadrature points
nQuad = size(lambda,1);


%% Compute geometric quantities and gradient of local basis
[Dphi,area] = gradbasis(node,elem);


%% ray information
ray_num = zeros(N,1);     % number of rays at each grid point
ray_dof = zeros(N,1);     % ray_dof(n) = sum(ray_num(1:n)) 

temp = 0;
for n = 1:N
    ray_num(n) = size(ray{n},2);
    ray_dof(n) = temp + ray_num(n);
    temp = ray_dof(n);
end

Ndof = temp;              % total degree of freedom

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


%% Assembling the linear system Av = b
rows = zeros(Nvec,1);
cols = rows;
vals = rows;
inds = 1;

b = zeros(Ndof,1);   

fprintf('Assembling time: \n');
tic;
for nt = 1:NT    
    for p = 1:nQuad
        % quadrature points in the x-y coordinate
        pxy = lambda(p,1)*node(elem(nt,1),:) ...
            + lambda(p,2)*node(elem(nt,2),:) ...
            + lambda(p,3)*node(elem(nt,3),:);
        fp = pde.f(pxy);
        k2 = (omega./pde.speed(pxy)).^2;
        
        sx = s_x(pxy(1),pxy(2));
        sy = s_y(pxy(1),pxy(2));
        sxy = s_xy(pxy(1),pxy(2));
        
        for i = 1:3
            ei = elem(nt,i);    % the ei'th grid point corresponding to the i'th vertix of element nt
            ni = ray_num(ei);   % number of rays crossing at grid point ei
            for nii = 1: ni
                gradtempi = - ray{ei}(nii);
                gradtempi = [real(gradtempi), imag(gradtempi)];
                fphasei = gradtempi(1)*pxy(1) + gradtempi(2)*pxy(2);
                ki = k(ei);
                phasei = exp(1i*ki*fphasei);
                
                ii = ray_dof(ei) - ni + nii;
                b(ii) = b(ii) + weight(p)*phi(p,i)*phasei*fp*sxy*area(nt);   % right hand side
                
                for j = 1:3
                    ej = elem(nt,j);    % the ej'th grid point corresponding to the j'th vertix of element nt
                    nj = ray_num(ej);   % number of rays crossing at grid point ej
                    
                    for njj = 1: nj
                        gradtempj = ray{ej}(njj);
                        gradtempj = [real(gradtempj), imag(gradtempj)];
                        fphasej = gradtempj(1)*pxy(1) + gradtempj(2)*pxy(2);
                        kj = k(ej);
                        phasej = exp(1i*kj*fphasej);
                        exp_phase = phasei*phasej;
                        
                        jj = ray_dof(ej) - nj + njj;
                        
                        tempA1 = sx*Dphi(nt,1,i)*Dphi(nt,1,j) + sy*Dphi(nt,2,i)*Dphi(nt,2,j);
                        
                        tempA2 = 1i*ki*phi(p,i)*(sx*gradtempi(1)*Dphi(nt,1,j) + sy*gradtempi(2)*Dphi(nt,2,j))...
                            + 1i*kj*phi(p,j)*(sx*gradtempj(1)*Dphi(nt,1,i) + sy*gradtempj(2)*Dphi(nt,2,i));
                        
                        tempA3 = ki*kj*phi(p,i)*phi(p,j)...
                            *(sx*gradtempi(1)*gradtempj(1) + sy*gradtempi(2)*gradtempj(2));
                        
                        tempM = k2*sxy*phi(p,i)*phi(p,j);
                        
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
toc;


%% Boundary conditions
[bdnode,~,~] = findboundary(elem);
isFreeNode = ones(Ndof,1);
Nbd = length(bdnode);
for nn = 1:Nbd
    nbd = bdnode(nn);
    ind = ray_dof(nbd) - ray_num(nbd) + 1:1:ray_dof(nbd);
    isFreeNode(ind) = 0;
end
freeNode = find(isFreeNode);


%% Solve the linear system Au = b
v = zeros(Ndof,1);
fprintf('Direct solver time:\n');
tic;
v(freeNode) = A(freeNode,freeNode)\b(freeNode);
toc;


%% Compute the solution value at grid nodes
u = zeros(N,1);
for n = 1:N
    grad = ray{n};
    grad = transpose(grad);
    grad = [real(grad), imag(grad)];    
    tempn = node(n,1)*grad(:,1) + node(n,2)*grad(:,2);
    ni = ray_dof(n) - ray_num(n);
    nii = 1:ray_num(n);
    nii = nii' + ni;
    u(n) = sum(v(nii).*exp(1i*k(n)*tempn));
end


%% Plot the solution
if plt    
    figure(16);
    FJ_showresult(node,elem,real(u));
    title('Ray-based FEM solution with PML');
end

