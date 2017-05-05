%% test for three interpolation methods:

a = 3;  h_c = 1/100;  h = 1/400;

%% 1) interp2
[X,Y] = meshgrid(-a:h_c:a, -2*a:h_c:2*a);
% [X,Y] = meshgrid(-a:h_c:a);
V = sin(X).*cos(Y);

figure(1)
surf(X,Y,V)
title('Original Sampling');

[Xq,Yq] = meshgrid(-a:h:a, -2*a:h:2*a); 
% [Xq,Yq] = meshgrid(-a:h:a);
Vq = interp2(X,Y,V,Xq,Yq);

figure(2)
surf(Xq,Yq,Vq);
title('Linear Interpolation Using Finer Grid');

% %% 2) interpolation
% [node,elem] = squaremesh([-a,a,-a,a],h_c);
% u = sin(node(:,1)).*cos(node(:,2));
% [xy,~] = squaremesh([-a,a,-a,a],h);
% uint = interpolation(node,elem,xy,u);
% 
% d1 = Vq(:) - uint;
% norm(d1, inf )


%% 3) interpolation2
[node,elem] = squaremesh([-a,a,-2*a,2*a],h_c);
u = sin(node(:,1)).*cos(node(:,2));
xmin = min(node(:,1));    xmax = max(node(:,1));
ymin = min(node(:,2));    ymax = max(node(:,2));
x = xmin:h_c:xmax; y = ymin:h_c:ymax; 
[xy,~] = squaremesh([-a,a,-2*a,2*a],h);
uint2 = interpolation2(x, y, u, xy);
d2 = Vq(:) - uint2;
norm(d2, inf )

figure(3); surf(Xq,Yq,reshape(uint2, size(Xq)));

% d3 = uint - uint2;
% norm(d3, inf )

