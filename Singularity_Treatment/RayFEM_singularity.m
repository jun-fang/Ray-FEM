function [u,A,b,v] = RayFEM_singularity(node,elem,xs,ys,omega,epsilon,wpml,sigmaMax,ray,speed,sing_rhs,fquadorder)
%% Ray-FEM solution with singularity treatment

A = assemble_Helmholtz_matrix_RayFEM(node,elem,omega,wpml,sigmaMax,speed,ray,fquadorder);
b = assemble_RHS_with_sing_RayFEM(node,elem,xs,ys,omega,epsilon,wpml,sigmaMax,ray,speed,sing_rhs,fquadorder);



%% Boundary conditions
[~,~,isBdNode] = findboundary(elem);
freeNode = find(~isBdNode);
N = size(node,1);  Nray = size(ray,2); 
v = zeros(N,1);
v(freeNode) = A(freeNode,freeNode)\b(freeNode);

%% construct solution
grad = ray(:);
grad = [real(grad),imag(grad)];
repnode = repmat(node,Nray,1);
temp = grad(:,1).*repnode(:,1) + grad(:,2).*repnode(:,2);

k = omega./speed(node);           % wavenumber
kk = repmat(k,1,Nray);
u = v.*exp(1i*kk(:).*temp);
u = reshape(u,N,Nray);
u = sum(u,2);