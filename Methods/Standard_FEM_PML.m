function [u,A,b] = Standard_FEM_PML(node,elem,omega,wpml,sigmaMax,pde,fquadorder,solver,plt)
%% Standard FEM for solving Helmholtz equation with PML: 
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
%   wpml: Width of the PML
% 
%   sigmaMax: Maximun absorbtion
%
%   pde: Structure of functions containing the pde information  
%        like speed c(x), right hand side f(x)
% 
%   fquadorder: The order of numerical quadrature
%
%   solver: The solver of linear system Au = b
%          - 'DIR': Direct solver A\b
%          - 'HIF': Hierachical Interplotative Factorization
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
N = size(node,1);        % number of grid points
NT = size(elem,1);       % number of triangle elements
Ndof = N;                % degree of freedom 

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
phi = lambda;           % linear bases
nQuad = size(lambda,1);


%% Compute geometric quantities and gradient of local basis
[Dphi,area] = gradbasis(node,elem);


%% Assemble the stiffness matrix A
% rows = zeros(9*nQuad*NT,1);
% cols = rows;
% vals = rows;
% inds = 1: NT;

A = sparse(Ndof,Ndof);
bt = zeros(NT,3);       % the right hand side

for p = 1:nQuad
    % quadrature points in the x-y coordinate
    pxy = lambda(p,1)*node(elem(:,1),:) ...
        + lambda(p,2)*node(elem(:,2),:) ...
        + lambda(p,3)*node(elem(:,3),:);
    k2 = (omega./pde.speed(pxy)).^2;
    sx = s_x(pxy(:,1),pxy(:,2));
    sy = s_y(pxy(:,1),pxy(:,2));
    sxy = s_xy(pxy(:,1),pxy(:,2));
    fp = pde.f(pxy).*sxy;
%     fp = pde.f(pxy).*( pxy(:,1) < xmax - wpml ).*( pxy(:,1) > xmin + wpml )...
%         .*( pxy(:,2) < ymax - wpml ).*( pxy(:,2) > ymin + wpml );
    
    for i = 1:3
        bt(:,i) = bt(:,i) + weight(p)*phi(p,i)*fp;
        for j = 1:3
            % $Delta_{ij}|_{\tau} = \int_{\tau} ( s1/s2 \partial_1 \phi_i \partial_1 \phi_j + s2/s1 \partial_2 \phi_i \partial_2 \phi_j ) dxdy$           
            Deltaij = sx.*Dphi(:,1,i).*Dphi(:,1,j) + sy.*Dphi(:,2,i).*Dphi(:,2,j);
            
            % $M_{ij}|_{\tau} = \int_{\tau} (omega/c)^2/(s1*s2) \phi_i \phi_j  dxdy$
            Mij = k2.*sxy*phi(p,i)*phi(p,j); 
            
            Aij = weight(p)*(Deltaij - Mij).*area; 
            A = A + sparse(elem(:,i),elem(:,j),Aij,Ndof,Ndof);   %% save memory
            
%             rows(inds) = elem(:,i);
%             cols(inds) = elem(:,j);
%             vals(inds) = Aij;           
%             inds = inds + NT;                
        end
    end
end

% A = sparse(rows,cols,vals,Ndof,Ndof);

bt = bt.*repmat(area,1,3);
b = accumarray(elem(:),bt(:),[Ndof 1]);

clear Aij Deltaij Mij Dphi k2 fp pxy;
clear bt area;



%% Boundary conditions
[bdNode,~,isBdNode] = findboundary(elem);
freeNode = find(~isBdNode);


%% Solve linear system Au = b
u = zeros(N,1);

% Direct solver
if (strcmp(solver,'DIR'))    
    fprintf('Direct solver time:\n');
    tic;
    u(freeNode) = A(freeNode,freeNode)\b(freeNode);
    toc;
end

% HIF-DE
if (strcmp(solver,'HIF'))    
    bdidx = zeros(N,1);
    bdidx(bdNode) = 1;
    Tbd = spdiags(bdidx,0,N,N);
    T = spdiags(1-bdidx,0,N,N);
    AD = T*A*T + Tbd;
    
    b(bdNode) = u(bdNode);
    A = AD;
       
    occ = 8;
    rank_or_tol = 1e-6;
    skip = 2;
    symm = 'n';
    n = round(sqrt(N)) + 1;
    opts = struct('skip',skip,'symm',symm,'verb',0);
    
    fprintf('HIF-DE solver time:\n');
    tic;
    F = hifde2(A,n,occ,rank_or_tol,opts);    
    [u,~,~,~] = gmres(@(x)(A*x),b,[],1e-9,32,@(x)(hifde_sv(F,x)));
    toc;
end


%% plot the solution
if plt
    figure(12);
    FJ_showresult(node,elem,real(u));
    title('Standard FEM solution with PML');
end



%% Memory estimate    total memory size: 24N + 23NT ~ 70N    Bytes: 70N*8 = 560N bytes 
% A:         21N,   7N nonzeros   rows,colums,values
% b:         N
% bt:        3NT

% elem:      3NT
% node:      2N

% Dphi:      6NT
% area:      NT

% Aij:       NT
% Deltaij:   NT
% Mij:       NT

% pxy:       2NT
% fp:        NT
% k2:        NT
% sx:        NT
% sy:        NT
% sxy:       NT


%% non compressible terms
% A,b,elem,node,Dphi,Aij,pxy: 24N + 12NT ~ 48N
