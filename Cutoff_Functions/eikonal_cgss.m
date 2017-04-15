function [ray, T] = eikonal_cgss(S02, g0, node0, node)
%% The analytical solution of eikonal equation with 
% constant gradient of slowness squared(CGSS)
% 
%         S^2 = |\nabla T|^2 = S0^2 + 2*g0 \cdot (\x - \x0)
% 
% two branch of solutions: |S +/- S0| ( S^2 + S0^2 +/- SS0 ) / (3|g0|)
%                                          or \frac{2|\x-\x_0| ( S^2 + S0^2 +/- SS0 )}{3|S -/+ S0|}
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


%% branch one: T1 = |S + S0| ( S^2 + S0^2 + SS0 ) / (3|g0|)
T1 = (S.^3 + 2*S0*S2 + 2*S02*S + S0^3) / (3*norm(g0));
dT1dx = (3*S2 + 4*S0*S + 2*S02).*dSdx / (3*norm(g0));
dT1dy = (3*S2 + 4*S0*S + 2*S02).*dSdy / (3*norm(g0));


%% branch two: T2 = \frac{2|\x-\x_0| ( S^2 + S0^2 - SS0 )}{3 (S + S0)}
N = 2*r.*(S2 + S02 - S.*S0);
D = 3*( S + S0);  D2 = D.*D;
T2 = N./D;

dNdx = 2*(drdx.*(S2 + S02 - S.*S0) + r.*(2*S - S0).*dSdx);
dNdy = 2*(drdy.*(S2 + S02 - S.*S0) + r.*(2*S - S0).*dSdy);
dDdx = 3*dSdx;  dDdy = 3*dSdy;

dT2dx = (dNdx.*D - N.*dDdx)./ D2;
dT2dy = (dNdy.*D - N.*dDdy)./ D2;


%% two braches:
T = [T1, T2];
ray1 = atan2(dT1dy, dT1dx);  ray1 = exp(1i*ray1);
ray2 = atan2(dT2dy, dT2dx);  ray2 = exp(1i*ray2);
ray2( r < 10*eps) = 0;
ray = [ray1, ray2]; 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% positive sign
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% forward pass:
% sumS = S + S0;
% recpsum = S02./sumS;
% temp = S + recpsum;
% T = 2/3*r.*temp;    % T = 2/3*r.*( S + S02./( S + S0 ));
% 
% %% backpropagation/chain rule to compute the gradient/ray direction
% dTdr = 2/3*temp;    % dT/dr
% dTdtemp = 2/3*r;    % dT/dtemp
% dS = 1*dTdtemp;
% drecpsim = 1*dTdtemp;
% dsumS = -S02./(sumS.^2).*drecpsim;
% dS = dS + 1*dsumS;
% 
% dTdx = drdx.*dTdr;
% dTdy = drdy.*dTdr;
% 
% dTdx = dTdx + dSdx.*dS;
% dTdy = dTdy + dSdy.*dS;
% 
% dTdx(r < 10*eps) = 0;
% dTdy(r < 10*eps) = 0;


% %% gradient check
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




