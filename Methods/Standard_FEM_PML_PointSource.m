function [u,A,b,h] = Standard_FEM_PML_PointSource(node,elem,h,omega,wpml,sigmaMax,xs,ys,speed,fquadorder,plt)
%% Standard FEM with PML for point source Helmholtz equation: 
%         -\Delta u - (omega/c)^2 u = \delta(x_0)         in D
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
%   wpml: Width of the PML
% 
%   sigmaMax: Maximun absorbtion
%
%   xs,ys: point source location
%
%   speed: wavespeed
% 
%   fquadorder: The order of numerical quadrature
%
%   plt: 1 plot the solution; 0 not plot
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% OUTPUT:
%   
%   u: N x 1 vector, solution of the linear system Au = b
%   A: N x N assembling matrix 
%   b: N x 1 vector, the right hand side
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% FEM set up
N = size(node,1);           % number of grid points
NT = size(elem,1);          % number of triangle elements
Ndof = N;                   % degree of freedom 

xmin = node(1,1);    xmax = node(end,1);
ymin = node(1,2);    ymax = node(end,2);
% h = (xmax-xmin)/round(sqrt(N));    % mesh size


%% PML set up
sigmaPML_x = @(x) sigmaMax*( (x-xmin-wpml).^2.*(x < xmin + wpml) + ...
                (x-(xmax-wpml)).^2.*(x > xmax - wpml))/wpml^2;         
sigmaPML_y = @(y) sigmaMax*( (y-ymin-wpml).^2.*(y < ymin + wpml) ...
                + (y-(ymax-wpml)).^2.*(y > ymax - wpml))/wpml^2;
s_x = @(x,y) (1+1i*sigmaPML_y(y)/omega)./(1+1i*sigmaPML_x(x)/omega);       %% s1/s2
s_y = @(x,y) (1+1i*sigmaPML_x(x)/omega)./(1+1i*sigmaPML_y(y)/omega);       %% s2/s1
s_xy = @(x,y) ((1+1i*sigmaPML_x(x)/omega).*(1+1i*sigmaPML_y(y)/omega));    %% 1/(s1*s2)


%% Numerical Quadrature
[lambda,weight] = quadpts(fquadorder);
phi = lambda;           % linear bases
nQuad = size(lambda,1);


%% Compute geometric quantities and gradient of local basis
[Dphi,area] = gradbasis(node,elem);


%% Assemble the stiffness matrix A
A = sparse(Ndof,Ndof);
bt = zeros(NT,3);       % the right hand side

for p = 1:nQuad
    % quadrature points in the x-y coordinate
    pxy = lambda(p,1)*node(elem(:,1),:) ...
        + lambda(p,2)*node(elem(:,2),:) ...
        + lambda(p,3)*node(elem(:,3),:);
    k2 = (omega./speed(pxy)).^2;
    sx = s_x(pxy(:,1),pxy(:,2));
    sy = s_y(pxy(:,1),pxy(:,2));
    sxy = s_xy(pxy(:,1),pxy(:,2));
    
    fp = nodal_basis(xs,ys,h,pxy);
   
    for i = 1:3
        bt(:,i) = bt(:,i) + weight(p)*phi(p,i)*fp;
        for j = 1:3
            % $Delta_{ij}|_{\tau} = \int_{\tau} ( s1/s2 \partial_1 \phi_i \partial_1 \phi_j + s2/s1 \partial_2 \phi_i \partial_2 \phi_j ) dxdy$           
%             Deltaij = sx.*Dphi(:,1,i).*Dphi(:,1,j) + sy.*Dphi(:,2,i).*Dphi(:,2,j);
%             
%             % $M_{ij}|_{\tau} = \int_{\tau} (omega/c)^2/(s1*s2) \phi_i \phi_j  dxdy$
%             Mij = sxy*phi(p,i)*phi(p,j); 
            
%             Aij = weight(p)*(Deltaij - k2.*Mij).*area;
%             MMij = weight(p)*Mij.*area;
            
            Aij = weight(p)*(sx.*Dphi(:,1,i).*Dphi(:,1,j) + sy.*Dphi(:,2,i).*Dphi(:,2,j) - k2.*sxy*phi(p,i)*phi(p,j)).*area;
            A = A + sparse(elem(:,i),elem(:,j),Aij,Ndof,Ndof);   %% save memory
        end
    end
end

bt = bt.*repmat(area,1,3);
b = accumarray(elem(:),bt(:),[N 1]);
b = b/(h*h/2);

%% Boundary conditions
[~,~,isBdNode] = findboundary(elem);
freeNode = find(~isBdNode);


%% Solve linear system Au = b
u = zeros(N,1);
% Direct solver
u(freeNode) = A(freeNode,freeNode)\b(freeNode);

%% plot the solution
if plt
    figure(12);
    FJ_showresult(node,elem,real(u));
    title('Standard FEM solution with PML');
end