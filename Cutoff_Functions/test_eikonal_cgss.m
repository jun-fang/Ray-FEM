%% Test for the analytical solution of eikonal equation 
% with constan gradient squared slowness  eikonal_cgss.m 
% 
% compare the results with example 2 in S. Luo and J. Qian's paper:
% http://users.math.msu.edu/users/qian/papers/LuoQianJCP2011.pdf

S02 = 4; g0 = [0, -3]; node0 = [0.25, 0];
h = 1/400;
[node, ~] = squaremesh([0, 0.5, -0.25, 0.5], h);
[~, T ] = eikonal_cgss(S02, g0, node0, node);

% figure(11);
% showsolution(node, elem, T(:,2), 2); colorbar;
% 
% figure(12);
% showsolution(node, elem, dT2dx, 2); colorbar;
% 
% figure(13);
% showsolution(node, elem, dT2dy, 2); colorbar;

x = -0:h:0.5;
y = -0.25:h:0.5;
[X,Y] = meshgrid(x,y);

figure(14);
Z = reshape(T(:,2), size(X));
contour(X,Y,Z, 20); colorbar;

[node, elem] = squaremesh([0, 0.5, 0.025, 0.5], h);
[ray, T, dT1dx, dT1dy, dT2dx, dT2dy] = eikonal_cgss(S02, g0, node0, node);

x = -0:h:0.5;
y = 0.025:h:0.5;
[X,Y] = meshgrid(x,y);


figure(15);
Z = reshape(dT2dx, size(X));
contour(X,Y,Z, 50); colorbar;


figure(16);
Z = reshape(dT2dy, size(X));
contour(X,Y,Z, 50); colorbar;

