global omega a xs ys;
pde = Helmholtz_data2;
xs = 2; ys = 2;
NPW = 6;
a = 1/2;
k = 4;

% omegas = [10 20 40 80]*pi;
% errors = 0*omegas;
% tic
% for i = 1:length(omegas)
%     i
%     tic;
%     omega = omegas(i);
%     h = 1/(NPW*round(omega/(2*pi)));
%     [node,elem] = squaremesh([-a,a,-a,a],h);
%     ray = pde.ray(node);
%     [error] = best_approximation(node,elem,ray,omega,pde);
%     errors(i) = error;
%     toc;
% end
% toc;
% 
% showrate(omegas,errors)

omega = 10*pi;
h = 1/(NPW*round(omega/(2*pi)));
[lnode,lelem] = squaremesh([-a,a,-a,a],h);
node = lnode(lelem(1,:),:);
lray = pde.ray(lnode);
ray = lray(lelem(1,:),:);

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

X1 = X(:,[1,5,9]);
y1 = pde.ex_u1(node);
c1 = transpose(X1)*X1 \ (transpose(X1)*y1);
err1 = norm(y1 - X1*c1)*h;

X2 = X(:,1+[1,5,9]);
y2 = pde.ex_u2(node);
c2 = transpose(X2)*X2 \ (transpose(X2)*y2);
err2 = norm(y2 - X2*c2)*h;

X3 = X(:,2+[1,5,9]);
y3 = pde.ex_u3(node);
c3 = transpose(X3)*X3 \ (transpose(X3)*y3);
err3 = norm(y3 - X3*c3)*h;

X4 = X(:,3+[1,5,9]);
y4 = pde.ex_u4(node);
c4 = transpose(X4)*X4 \ (transpose(X4)*y4);
err4 = norm(y4 - X4*c4)*h;



c11 = pde.ex_u1(node1)/exp(1i*omega* ( real(ray(1,1))*node1(1) + imag(ray(1,1))*node1(2) ));
c12 = pde.ex_u1(node2)/exp(1i*omega* ( real(ray(2,1))*node2(1) + imag(ray(2,1))*node2(2) ));
c13 = pde.ex_u1(node3)/exp(1i*omega* ( real(ray(3,1))*node3(1) + imag(ray(3,1))*node3(2) ));

cc1 = transpose([c11,c12,c13]);
err11 = norm(y1 - X1*cc1)*h;





