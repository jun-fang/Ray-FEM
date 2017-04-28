%% Test for the analytical solution of eikonal equation 
% with constan gradient squared slowness: eikonal_cgss.m 
% and constant gradient of velocity: eikonal_cgv.m
% 
% compare the results with example 2 and 3 in S. Luo and J. Qian's paper:
% http://users.math.msu.edu/users/qian/papers/LuoQianJCP2011.pdf


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 2: constan gradient squared slowness: eikonal_cgss.m 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% gradient check
% S02 = 1;
% g0 = [1, 2];
% node0 = [0.32, -0.0851];
% h = 1e-6;
% node = [0.5,0.5];
% [ray, T, dT1dx, dT1dy, dT2dx, dT2dy] = eikonal_cgss(S02, g0, node0, node);
% 
% node1 = [0.5+h,0.5];
% node2 = [0.5-h,0.5];
% [~, T1] = eikonal_cgss(S02, g0, node0, node1);
% [~, T2] = eikonal_cgss(S02, g0, node0, node2);
% 
% (T1(:,1) - T2(:,1))/(2*h) - dT1dx
% (T1(:,2) - T2(:,2))/(2*h) - dT2dx
% 
% 
% node3 = [0.5,0.5+h];
% node4 = [0.5,0.5-h];
% [~, T3] = eikonal_cgss(S02, g0, node0, node3);
% [~,T4] = eikonal_cgss(S02, g0, node0, node4);
% 
% (T3(:,1) - T4(:,1))/(2*h) - dT1dy
% (T3(:,2) - T4(:,2))/(2*h) - dT2dy
% 


%% jianliang's paper
% S02 = 4; g0 = [0, -3]; node0 = [0.25, 0];
% h = 1/400;
% [node, ~] = squaremesh([0, 0.5, -0.25, 0.5], h);
% [~, T ] = eikonal_cgss(S02, g0, node0, node);
% 
% % figure(11);
% % showsolution(node, elem, T(:,2), 2); colorbar;
% % 
% % figure(12);
% % showsolution(node, elem, dT2dx, 2); colorbar;
% % 
% % figure(13);
% % showsolution(node, elem, dT2dy, 2); colorbar;
% 
% x = -0:h:0.5;
% y = -0.25:h:0.5;
% [X,Y] = meshgrid(x,y);
% 
% figure(14);
% Z = reshape(T(:,2), size(X));
% contour(X,Y,Z, 20); colorbar;
% 
% [node, elem] = squaremesh([0, 0.5, 0.025, 0.5], h);
% [ray, T, dT1dx, dT1dy, dT2dx, dT2dy] = eikonal_cgss(S02, g0, node0, node);
% 
% x = -0:h:0.5;
% y = 0.025:h:0.5;
% [X,Y] = meshgrid(x,y);
% 
% 
% figure(15);
% Z = reshape(dT2dx, size(X));
% contour(X,Y,Z, 50); colorbar;
% 
% 
% figure(16);
% Z = reshape(dT2dy, size(X));
% contour(X,Y,Z, 50); colorbar;



%% verify the eikonal equation
% % S02 = 4; g0 = [0, 1/1]; 
% node0 = [0, 0];
% omega0 = 40*pi; 0.8;
% E0 = omega0^2;                   % E0 should be large enough to make sure the medium is smooth enough
% speed = @(p) omega0./sqrt( E0 + p(:,2) );    % wave speed
% 
% 
% S02 = ( 1/speed(node0) )^2;
% g0 = [0, 1/(2*omega0*omega0)];
% % g0 = [0.1, 0.2];
% 
% h = 1/80;
% 
% [node, elem] = squaremesh([-0.5, 0.5, -0.5, 0.5], h);
% [ray, T, dT1dx, dT1dy, dT2dx, dT2dy] = eikonal_cgss(S02, g0, node0, node);
% 
% figure(17); 
% ray_field(ray,node,10,1/10);
% 
% x = node(:,1); y = node(:,2);
% x0 = node0(1,1); y0 = node0(1,2);
% S2 = S02 + 2* ( g0(1)*(x-x0) + g0(2)*(y-y0) );
% 
% dT1 = dT1dx.^2 +  dT1dy.^2;   df1 = dT1 - S2;
% dT2 = dT2dx.^2 +  dT2dy.^2;   df2 = dT2 - S2;
% sum (dT1 - S2)
% sum (dT2 - S2)



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 3: constant gradient of velocity: eikonal_cgv.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% gradient check
v0 = 2;
G0 = [0, -1];
node0 = [0, 0];
h = 1e-6;
x = -0.2; y = -0.3;
node = [x,y];
[ray, T, dTdx, dTdy] = eikonal_cgv(v0, G0, node0, node);

node1 = [x+h,y];
node2 = [x-h,y];
[~, T1] = eikonal_cgv(v0, G0, node0, node1);
[~, T2] = eikonal_cgv(v0, G0, node0, node2);

(T1 - T2)/(2*h) - dTdx


node3 = [x,y+h];
node4 = [x,y-h];
[~, T3] = eikonal_cgv(v0, G0, node0, node3);
[~, T4] = eikonal_cgv(v0, G0, node0, node4);

(T3 - T4)/(2*h) - dTdy



%% verify the eikonal equation
% S02 = 4; g0 = [0, 1/1]; 

v0 = 1;
G0 = [0, -1/10];
node0 = [0, 0];

h = 1/180;
[node, elem] = squaremesh([-0.5, 0.5, -0.5, 0.5], h);
[ray, T, dTdx, dTdy] = eikonal_cgv(v0, G0, node0, node);

figure(17); 
ray_field(ray,node,10,1/10);

x = node(:,1); y = node(:,2);
x0 = node0(1,1); y0 = node0(1,2);
v = v0 + ( G0(1)*(x-x0) + G0(2)*(y-y0) );

dT2 = dTdx.^2 +  dTdy.^2;   df = dT2 - 1./(v.*v);
r = sqrt( (x-x0).^2 + (y-y0).^2 );       % |\x - \x0|
df( r < 10*eps) = 0;
norm (df, inf)


