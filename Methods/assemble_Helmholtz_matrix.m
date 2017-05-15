function A = assemble_Helmholtz_matrix(node,elem,omega,wpml,sigmaMax,speed,fquadorder)
%% Function to assemble the standard FEM Helmholtz matrix with PML: 
%         -\Delta u - (omega/c)^2 u = f               in D
%                                 u = 0               on \partial D 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%   wpml: vector (4,1) Width of the PML in each direction 
% 
%   sigmaMax: Maximun absorbtion
%
%   pde: Structure of functions containing the pde information  
%        like speed c(x), right hand side f(x)
% 
%   fquadorder: The order of numerical quadrature
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% OUTPUT:
%   
%   A: N x N Helmholtx matrix
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% FEM set up
fprintf('Assembling Standard Finite Element matrix \n')

N = size(node,1);        % number of grid points
%NT = size(elem,1);       % number of triangle elements
Ndof = N;                % degree of freedom 

% computing the limits ofthe domain
xmax = max(node(:,1));
xmin = min(node(:,1));
ymax = max(node(:,2));
ymin = min(node(:,2));

% in order to define the correct boundary conditionts at the interfaces
% we need to specify the width of the PML in eevry direction
%    ---wpml(1)---
%   |             |
%   |             |
%  wpml(3)     wpml(4)
%   |             |
%   |             |
%    ---wpml(2)---
%

% if the different lengths of the pml are not specified we specify all of
% them uisng the same number of pml points

if length(wpml) == 1
    wpml = ones(4,1)*wpml;
end

%% PML set up

% % usual quadratic profile 
% sigmaPML_x = @(x)sigmaMax*( (x-xmin-wpml(1)).^2.*(x < xmin + wpml(1))/wpml(1)^2 + ...
%                 (x-(xmax-wpml(2))).^2.*(x > xmax - wpml(2))/wpml(2)^2);         
% sigmaPML_y = @(y) sigmaMax*( (y-ymin-wpml(3)).^2.*(y < ymin + wpml(3))/wpml(3)^2 ...
%                 + (y-(ymax-wpml(4))).^2.*(y > ymax - wpml(4))/wpml(4)^2);

% unbounded PML profile             
sigmaPML_x = @(p) speed(p).*(1./(p(:,1)-xmin + sigmaMax ).*(p(:,1) < xmin + wpml(1)) + ...
                1./(xmax-p(:,1)  + sigmaMax ).*(p(:,1) > xmax - wpml(2))) ;         
sigmaPML_y = @(p) speed(p).*( 1./(p(:,2)-ymin  + sigmaMax).*(p(:,2) < ymin + wpml(3)) ...
                + 1./(ymax-p(:,2) + sigmaMax).*(p(:,2) > ymax - wpml(4)));
             
            
s_x = @(p) (1+1i*sigmaPML_y(p)/omega)./(1+1i*sigmaPML_x(p)/omega);       %% s1/s2
s_y = @(p) (1+1i*sigmaPML_x(p)/omega)./(1+1i*sigmaPML_y(p)/omega);       %% s2/s1
s_xy = @(p) ((1+1i*sigmaPML_x(p)/omega).*(1+1i*sigmaPML_y(p)/omega));    %% 1/(s1*s2)


%% Numerical Quadrature
[lambda,weight] = quadpts(fquadorder);
phi = lambda;           % linear bases
nQuad = size(lambda,1);


%% Compute geometric quantities and gradient of local basis
[Dphi,area] = gradbasis(node,elem);


%% Assembling the matrix

A = sparse(Ndof,Ndof);

for p = 1:nQuad
    % quadrature points in the x-y coordinate
    pxy =   lambda(p,1)*node(elem(:,1),:) ...
          + lambda(p,2)*node(elem(:,2),:) ...
          + lambda(p,3)*node(elem(:,3),:);
      
    % local wavenumber
    k2  = (omega./speed(pxy)).^2;
    % building the PML profiles
    sx  = s_x(pxy);
    sy  = s_y(pxy);
    sxy = s_xy(pxy);
    
    for i = 1:3
        for j = 1:3
            % $Delta_{ij}|_{\tau} = \int_{\tau} ( s1/s2 \partial_1 \phi_i \partial_1 \phi_j + s2/s1 \partial_2 \phi_i \partial_2 \phi_j ) dxdy$           
            Deltaij = sx.*Dphi(:,1,i).*Dphi(:,1,j) + sy.*Dphi(:,2,i).*Dphi(:,2,j);
            
            % $M_{ij}|_{\tau} = \int_{\tau} (omega/c)^2/(s1*s2) \phi_i \phi_j  dxdy$
            Mij = k2.*sxy*phi(p,i)*phi(p,j); 
            
            Aij = weight(p)*(Deltaij - Mij).*area; 
            A = A + sparse(elem(:,i),elem(:,j),Aij,Ndof,Ndof);   %% save memory
                          
        end
    end
end

clear Aij Deltaij Mij Dphi k2 fp pxy;
clear area;