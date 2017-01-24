function [dTdx, dTdy, T] = eikonal_cgss(S02, g0, node0, node)
%% The analytical solution of eikonal equation with 
% constant gradient of slowness squared(CGSS)
% 
%         |\nabla T|^2 = S0^2 + 2*g0 \cdot (\x - \x0)
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

%% forward pass:
sumS = S + S0;
recpsum = S02./sumS;
temp = S + recpsum;
T = 2/3*r.*temp;    % T = 2/3*r.*( S + S02./( S + S0 ));

%% backpropagation/chain rule to compute the gradient/ray direction
dTdr = 2/3*temp;    % dT/dr
dTdtemp = 2/3*r;    % dT/dtemp
dS = 1*dTdtemp;
drecpsim = 1*dTdtemp;
dsumS = -S02./(sumS.^2).*drecpsim;
dS = dS + 1*dsumS;

dTdx = drdx.*dTdr;
dTdy = drdy.*dTdr;

dTdx = dTdx + dSdx.*dS;
dTdy = dTdy + dSdy.*dS;

dTdx(r < 10*eps) = 0;
dTdy(r < 10*eps) = 0;


%% gradient check
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




