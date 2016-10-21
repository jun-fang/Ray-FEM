function b = assemble_ref_RHS(node,elem,h,xc,yc,fquadorder)
N = size(node,1);        % number of grid points
NT = size(elem,1);       % number of triangle elements

%% Numerical Quadrature
[lambda,weight] = quadpts(fquadorder);
phi = lambda;           % linear bases
nQuad = size(lambda,1);


%% Compute geometric quantities and gradient of local basis
[~,area] = gradbasis(node,elem);

%% Assemble right-hand side

bt = zeros(NT,3);       % the right hand side

for p = 1:nQuad
    % quadrature points in the x-y coordinate
    pxy = lambda(p,1)*node(elem(:,1),:) ...
        + lambda(p,2)*node(elem(:,2),:) ...
        + lambda(p,3)*node(elem(:,3),:);
    
    % we suppose that the source is well inside the physical domain
    fp = approximate_delta(xc,yc,h,pxy);
%     fp = source(pxy).*( pxy(:,1) < xmax - wpml ).*( pxy(:,1) > xmin + wpml )...
%         .*( pxy(:,2) < ymax - wpml ).*( pxy(:,2) > ymin + wpml ); 
    for i = 1:3
        bt(:,i) = bt(:,i) + weight(p)*phi(p,i)*fp;
    end
end

bt = bt.*repmat(area,1,3);
b = accumarray(elem(:),bt(:),[N 1]);
