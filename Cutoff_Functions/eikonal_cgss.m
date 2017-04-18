function [ray, T, dT1dx, dT1dy, dT2dx, dT2dy] = eikonal_cgss(S02, g0, node0, node)
%% The analytical solution of eikonal equation with 
% constant gradient of slowness squared(CGSS)
% 
%         S^2 = |\nabla T|^2 = S0^2 + 2*g0 \cdot (\x - \x0)
% 
% two branch of solutions: |S +/- S0| ( S^2 + S0^2 -/+ SS0 ) / (3|g0|)
%                                          or \frac{2|\x-\x_0| ( S^2 + S0^2 -/+ SS0 )}{3|S -/+ S0|}
% 
% S02: 1x1
% g0: 1x2
% node0: 1x2
% node: Nx2

x = node(:,1);  y = node(:,2);
x0 = node0(:,1);  y0 = node0(:,2);
r = sqrt((x-x0).^2 + (y-y0).^2);       % |\x - \x0|
drdx = (x-x0)./r;  drdy = (y-y0)./r;   % dr/dx, dr/dy 

S2 = S02 + 2*( g0(1)*(x-x0) + g0(2)*(y-y0));
S = sqrt(S2); S0 = sqrt(S02);
dSdx = g0(1)./S;  dSdy = g0(2)./S;


%% branch one: T1 = (S^3 + s0^3) / (3|g0|)
gr = sqrt(g0(1)^2 + g0(2)^2);
T1 = (S.^3 + S0^3) / (3*gr);
dT1dx = S2.*dSdx / gr;
dT1dy = S2.*dSdy / gr;


%% branch two: T2 =2/3 * |\x-\x_0| ( S + S0^2 / (S + S0) )
F =  S + S02./ (S + S0);
T2 = 2/3*r.*F;

dF = 1 - S02./ (S+S0).^2;
dT2dx = 2/3*( drdx.*F + r.*dF.*dSdx );
dT2dy = 2/3*( drdy.*F + r.*dF.*dSdy );


%% two braches:
T = [T1, T2];
ray1 = atan2(dT1dy, dT1dx);  ray1 = exp(1i*ray1);
ray2 = atan2(dT2dy, dT2dx);  ray2 = exp(1i*ray2);
ray2( r < 10*eps) = 0;
ray = [ray1, ray2]; 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gradient check
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  
% S02 = 1;
% g0 = [1, 2];
% node0 = [-0.32, -0.111];
% h = 1e-6;
% node = [0.5,0.5];
% [dx, dy, T] = eikonal_cgss(S02, g0, node0, node);
% 
% node1 = [0.5+h,0.5];
% node2 = [0.5-h,0.5];
% [~,~, T1] = eikonal_cgss(S02, g0, node0, node1);
% [~,~, T2] = eikonal_cgss(S02, g0, node0, node2);
% 
% (T1 - T2)/(2*h) - dx
% 
% node3 = [0.5,0.5+h];
% node4 = [0.5,0.5-h];
% [~,~, T3] = eikonal_cgss(S02, g0, node0, node3);
% [~,~, T4] = eikonal_cgss(S02, g0, node0, node4);
% 
% (T3 - T4)/(2*h) - dy




