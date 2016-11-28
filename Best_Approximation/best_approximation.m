function [error] = best_approximation(node,elem,ray,omega,pde)
%% compute a lower bound of the best approximation in the ray-FEM space
NT = size(elem,1);
error = 0;
for i = 1:NT
    node_i = node(elem(i,:),:);
    ray_i = ray(elem(i,:),:);
    [~,err] = best_approximation_in_one_element(node_i,ray_i,omega,pde,4);
    error = error + err*err;
end
error = sqrt(error);
end


function [c, err] = best_approximation_in_one_element(node,ray,omega,pde,k)
%% optimize the approximation error in one element by least squares

elem = [1 2 3];
node1 = node(1,:); node2 = node(2,:); node3 = node(3,:);
h = max(abs(node1-node2));
elem_area = tri_area(node1, node2, node3);

for i = 1:k
    [node,elem] = uniformrefine(node,elem);
    h = h/2;
end
N = size(node,1); Nray = size(ray,2);

phi1 = tri_area(node, repmat(node2,N,1), repmat(node3,N,1))/elem_area;
phi2 = tri_area(node, repmat(node3,N,1), repmat(node1,N,1))/elem_area;
phi3 = tri_area(node, repmat(node1,N,1), repmat(node2,N,1))/elem_area;
phi1 = repmat(phi1,1,Nray);
phi2 = repmat(phi2,1,Nray);
phi3 = repmat(phi3,1,Nray);

phase1 = zeros(N,Nray); phase2 = zeros(N,Nray); phase3 = zeros(N,Nray);
for ni = 1:Nray
    phase1(:,ni) = real(ray(1,ni)).*node(:,1) + imag(ray(1,ni)).*node(:,2);
    phase2(:,ni) = real(ray(2,ni)).*node(:,1) + imag(ray(2,ni)).*node(:,2);
    phase3(:,ni) = real(ray(3,ni)).*node(:,1) + imag(ray(3,ni)).*node(:,2);
end

X = [phi1.*exp(1i*omega*phase1), phi2.*exp(1i*omega*phase2), phi3.*exp(1i*omega*phase3)];
y = pde.ex_u(node);
% c = transpose(X)*X \ (transpose(X)*y);
% err = norm(y - X*c)*h;

% A = [real(X), -imag(X); imag(X), real(X)];
% b = [real(y); imag(y)];
% x = (A'*A)\(A'*b);
% nx = length(x);
% cc = x(1:nx/2) + 1i*x(nx/2+1:end);
% err = norm(abs(y - X*cc))*h;
% c = cc;

c11 = pde.ex_u1(node1)/exp(1i*omega* ( real(ray(1,1))*node1(1) + imag(ray(1,1))*node1(2) ));
c12 = pde.ex_u2(node1)/exp(1i*omega* ( real(ray(1,2))*node1(1) + imag(ray(1,2))*node1(2) ));
c13 = pde.ex_u3(node1)/exp(1i*omega* ( real(ray(1,3))*node1(1) + imag(ray(1,3))*node1(2) ));
c14 = pde.ex_u4(node1)/exp(1i*omega* ( real(ray(1,4))*node1(1) + imag(ray(1,4))*node1(2) ));

c21 = pde.ex_u1(node2)/exp(1i*omega* ( real(ray(2,1))*node2(1) + imag(ray(2,1))*node2(2) ));
c22 = pde.ex_u2(node2)/exp(1i*omega* ( real(ray(2,2))*node2(1) + imag(ray(2,2))*node2(2) ));
c23 = pde.ex_u3(node2)/exp(1i*omega* ( real(ray(2,3))*node2(1) + imag(ray(2,3))*node2(2) ));
c24 = pde.ex_u4(node2)/exp(1i*omega* ( real(ray(2,4))*node2(1) + imag(ray(2,4))*node2(2) ));

c31 = pde.ex_u1(node3)/exp(1i*omega* ( real(ray(3,1))*node3(1) + imag(ray(3,1))*node3(2) ));
c32 = pde.ex_u2(node3)/exp(1i*omega* ( real(ray(3,2))*node3(1) + imag(ray(3,2))*node3(2) ));
c33 = pde.ex_u3(node3)/exp(1i*omega* ( real(ray(3,3))*node3(1) + imag(ray(3,3))*node3(2) ));
c34 = pde.ex_u4(node3)/exp(1i*omega* ( real(ray(3,4))*node3(1) + imag(ray(3,4))*node3(2) ));

c = [c11 c12 c13 c14 c21 c22 c23 c24 c31 c32 c33 c34]';
err = norm(abs(y - X*c))*h;


end


