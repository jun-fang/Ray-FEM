function T = eikonal_cgss(S0, grad0, node0, node)
%% The analytical solution of eikonal equation with 
% constant gradient of slowness squared(CGSS)
% 
%         |\nabla T|^2 = S0^2 + 2*grad0 \cdot (\x - \x0)
% 
x = node(:,1);  y = node(:,2);
x0 = node0(:,1);  y0 = node0(:,2);
r = sqrt((x-x0).^2 + (y-y0).^2);

S2 = S0^2 + 2*( grad0(1)*(x-x0) + grad0(2)*(y-y0));
S = sqrt(S2);

T = 2/3*r.*( S2 + S0^2 + S.*S0 )./( S + S0 );