function [u,A,v,b] = Ray_FEM_PML_1_PointSource(node,elem,omega,wpml,sigmaMax,xs,ys,speed,ray,fquadorder,plt)
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
Nray = size(ray,2);     % number of rays crossing at each grid node
Ndof = N*Nray;          % degree of freedom 
% n = round(sqrt(N));         % number of grid points in each dimension

c = speed(node);    % medium speed
k = omega./c;           % wavenumber
kk = repmat(k,1,Nray);

xmin = node(1,1);    xmax = node(end,1);
ymin = node(1,2);    ymax = node(end,2);

h = node(2,2) - node(1,2);    % mesh size
n = round((ymax-ymin)/h + 1);

xn = round((xs - xmin)/h);
yn = round((ys - ymin)/h);
xyn = n*xn + yn + 1;


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


%% Assemble the linear system Av = b
A = sparse(Ndof,Ndof);
M = A;

fprintf('Assembling time:\n');
tic;
for p = 1:nQuad
    pxy = lambda(p,1)*node(elem(:,1),:) ...
        + lambda(p,2)*node(elem(:,2),:) ...
        + lambda(p,3)*node(elem(:,3),:);
    sx = s_x(pxy(:,1),pxy(:,2));
    sy = s_y(pxy(:,1),pxy(:,2));
    sxy = s_xy(pxy(:,1),pxy(:,2));
    k2 = (omega./speed(pxy)).^2;
    
    for i = 1:3       
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
                    
                    ii = (nii-1)*N + elem(:,i);
                    jj = (njj-1)*N + elem(:,j);
                    Aij = weight(p)*(tempA1 + tempA2 - tempA3 - tempM).*exp_phase.*area;
                    Mij = weight(p)*sxy*phi(p,i)*phi(p,j).*exp_phase.*area;
                    
                    A = A + sparse(ii,jj,Aij,Ndof,Ndof); 
                    M = M + sparse(ii,jj,Mij,Ndof,Ndof); 
                    
                end
            end
        end
    end
end
toc;

b = M(:,xyn)/h^2;
clear Aij Mij Dphi fphasei fphasej phasei phasej exp_phase gradtempi gradtempj;
clear k2 ki kj kki pxy fp reppxy tempA1 tempA2 tempA3 tempM;


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



