function [u,A,b,rel_L2_err] = Standard_FEM_IBC(node,elem,omega,pde,fquadorder,solver,plt)
%% Standard FEM for solving Helmholtz equation with Impedance Boundary Condition(IBC): 
%         -\Delta u - (omega/c)^2 u = f               in D
%         \partial u / \partial n + i omega/c u = g   on \partial D 
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
%   pde: Structure of functions containing the pde information  
%        like speed c(x), right hand side f(x)
% 
%   fquadorder: The order of numerical quadrature
%
%   solver: The solver of linear system Au = b
%          - 'Dir': Direct solver A\b
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
%   rel_L2_err: relative L2 error
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% FEM set up
N = size(node,1);       % number of grid points
NT = size(elem,1);      % number of triangle elements
Ndof = N;               % degree of freedom 


%% Numerical Quadrature
[lambda,weight] = quadpts(fquadorder);
phi = lambda;           % linear bases
nQuad = size(lambda,1);


%% Compute geometric quantities and gradient of local basis
[Dphi,area] = gradbasis(node,elem);


%% Assemble the linear system Au = b
% rows = zeros(nQuad*3*3*NT,1);
% cols = rows;
% vals = rows;
% inds = 1: NT;
A = sparse(Ndof,Ndof);
bt = zeros(NT,3);       % the right hand side

% fprintf('Assembling time:\n');
% tic;
for p = 1:nQuad
    % quadrature points in the x-y coordinate
    pxy = lambda(p,1)*node(elem(:,1),:) ...
        + lambda(p,2)*node(elem(:,2),:) ...
        + lambda(p,3)*node(elem(:,3),:);
    fp = pde.f(pxy);
    k2 = (omega./pde.speed(pxy)).^2;
    
    for i = 1:3         
        bt(:,i) = bt(:,i) + weight(p)*phi(p,i)*fp;
        
        for j = 1:3     % can be optimized because of the symmetry of (i,j)
            % $Delta_{ij}|_{\tau} = \int_{\tau} \nabla \phi_i \dot \nabla \phi_j dxdy$           
            Deltaij = Dphi(:,1,i).*Dphi(:,1,j) + Dphi(:,2,i).*Dphi(:,2,j);
            
            % $M_{ij}|_{\tau} = \int_{\tau} (omega/c)^2 \phi_i \phi_j  dxdy$
            Mij = k2.*phi(p,i)*phi(p,j);  
            Aij = weight(p)*(Deltaij - Mij).*area; 
            
            A = A + sparse(elem(:,i),elem(:,j),Aij,Ndof,Ndof);   %% save memory
            
%             rows(inds) = elem(:,i);        % maybe faster but cost much more memories
%             cols(inds) = elem(:,j);
%             vals(inds) = Aij;           
%             inds = inds + NT;                
        end
    end
end

% A = sparse(rows,cols,vals,Ndof,Ndof);
% clear Aij Deltaij Mij Dphi rows cols vals k2 fp pxy;
clear Aij Deltaij Mij Dphi k2 fp pxy;

bt = bt.*repmat(area,1,3);
b = accumarray(elem(:),bt(:),[Ndof 1]);

clear bt area;
% toc;


%% Boundaries
[~,bdEdge,~] = findboundary(elem);
el = sqrt(sum((node(bdEdge(:,1),:) - node(bdEdge(:,2),:)).^2,2));
Ne = length(el);

[lambdagN,weightgN] = quadpts1(fquadorder);
phigN = lambdagN;               
nQuadgN = size(lambdagN,1);


%% Modify the linear system Au = b with Impedance Boundary Condition
rows = zeros(nQuadgN*2*2*Ne, 1);
cols = rows;
vals = rows;
inds = 1:Ne;
ge = zeros(size(bdEdge,1),2);

for pp = 1:nQuadgN
    ppxy = lambdagN(pp,1)*node(bdEdge(:,1),:) ...
        + lambdagN(pp,2)*node(bdEdge(:,2),:);
    gp = pde.g_IBC(ppxy);
    kp = omega./pde.speed(ppxy);
    for i = 1:2
        ge(:,i) = ge(:,i) + weightgN(pp)*phigN(pp,i)*gp;
        for j = 1:2
            rows(inds) = bdEdge(:,i);
            cols(inds) = bdEdge(:,j);
            vals(inds) = weightgN(pp)*1i*phigN(pp,i)*phigN(pp,j)*kp.*el;
            inds = inds + Ne; 
        end
    end
end

A = A + sparse(rows,cols,vals,Ndof,Ndof);

ge = ge.*repmat(el,1,2);
b = b + accumarray(bdEdge(:), ge(:),[Ndof,1]);

clear rows cols vals; 

%% Solve the linear system Au = b
% Direct solver
if (strcmp(solver,'DIR'))    
%     fprintf('Direct solver time:\n');
%     tic;
    u = A\b;
%     toc;
end

% HIF-DE
if (strcmp(solver,'HIF'))    
    occ = 8;
    rank_or_tol = 1e-9;
    skip = 2;
    symm = 'n';
    n = round(sqrt(N)) + 1;
    opts = struct('skip',skip,'symm',symm,'verb',0);
    
    fprintf('HIF-DE solver time:\n');
    tic;
    F = hifde2(A,n,occ,rank_or_tol,opts);
    [u,~,~,~] = gmres(@(x)(A*x),b,[],1e-12,32,@(x)(hifde_sv(F,x)));
    toc;
end


%% Plot the solution
if plt
    figure(10);
    uex = pde.ex_u(node);
    FJ_showresult(node,elem,real(uex));
    title('Exact Solution');
    
    figure(11);
    FJ_showresult(node,elem,real(u));
    title('Standard FEM solution with IBC');
end


%% get relative L2 error
err = getL2error(node,elem,@pde.ex_u,u);
err0 = getL2error(node,elem,@pde.ex_u,0*u);
rel_L2_err = err/err0;

% uex = pde.ex_u(node);
% rel_L2_err = norm(u-uex)/norm(uex);



