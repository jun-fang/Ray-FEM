function [u, v] = RayFEM_direct_solver(node,elem,A,b,omega,ray,speed)
%% Direct solver of linear system Au = b with PML: u = A\b
% ray could be a N x Nray matrix or N x 1 cell

% Boundaries
[bdnode,~,isBdNode] = findboundary(elem);
N = size(node,1);  
k = omega./speed(node);             % wavenumber


%% ray is N x Nray matrix: each grid point has the same number of rays
if ~iscell(ray)
    
    freeNode = find(~isBdNode);
    Nray = size(ray,2);
    freeNode = repmat(freeNode, Nray, 1);
    v = zeros(N*Nray,1);
    v(freeNode) = A(freeNode,freeNode)\b(freeNode);
    
    % construct solution
    grad = ray(:);
    grad = [real(grad),imag(grad)];
    repnode = repmat(node,Nray,1);
    temp = grad(:,1).*repnode(:,1) + grad(:,2).*repnode(:,2);
    
    kk = repmat(k,Nray,1);
    u = v.*exp(1i*kk(:).*temp);
    u = reshape(u,N,Nray);
    u = sum(u,2);
    

else
    %% ray is N x 1 cell: each grid point may have different number of rays
    
    ray_num = zeros(N,1);     % number of rays at each grid point
    ray_dof = zeros(N,1);     % ray_dof(n) = sum(ray_num(1:n))
    
    temp = 0;
    for n = 1:N
        ray_num(n) = size(ray{n},2);
        ray_dof(n) = temp + ray_num(n);
        temp = ray_dof(n);
    end
    
    Ndof = temp;              % total degree of freedom
    
    isFreeNode = ones(Ndof,1);
    Nbd = length(bdnode);
    for nn = 1:Nbd
        nbd = bdnode(nn);
        ind = ray_dof(nbd) - ray_num(nbd) + 1:1:ray_dof(nbd);
        isFreeNode(ind) = 0;
    end
    freeNode = find(isFreeNode);
    
    % solve the linear system Au = b
    v = zeros(Ndof,1);
    v(freeNode) = A(freeNode,freeNode)\b(freeNode);    
 
    % construct solution
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
    
end
