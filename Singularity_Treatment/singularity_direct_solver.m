function u = singularity_direct_solver(node,elem,A,b,omega,ray,speed)
%% Direct solver of linear system Au = b with PML: u = A\b

% Boundary conditions
[~,~,isBdNode] = findboundary(elem);
freeNode = find(~isBdNode);
N = size(node,1);  Nray = size(ray,2); 
v = zeros(N,1);
v(freeNode) = A(freeNode,freeNode)\b(freeNode);

% construct solution
grad = ray(:);
grad = [real(grad),imag(grad)];
repnode = repmat(node,Nray,1);
temp = grad(:,1).*repnode(:,1) + grad(:,2).*repnode(:,2);

k = omega./speed(node);           % wavenumber
kk = repmat(k,1,Nray);
u = v.*exp(1i*kk(:).*temp);
u = reshape(u,N,Nray);
u = sum(u,2);