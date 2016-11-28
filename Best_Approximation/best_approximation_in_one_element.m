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
c = transpose(X)*X \ (transpose(X)*y);
err = norm(y - X*c)*h;

end

% c1 = pde.ex_u(node1)/exp(1i*omega* ( real(ray(1))*node1(1) + imag(ray(1))*node1(2) ));
% c2 = pde.ex_u(node2)/exp(1i*omega* ( real(ray(2))*node2(1) + imag(ray(2))*node2(2) ));
% c3 = pde.ex_u(node3)/exp(1i*omega* ( real(ray(3))*node3(1) + imag(ray(3))*node3(2) ));
% 
% c2 = [c1;c2;c3];
% err2 = norm(y - X*c2)*h;