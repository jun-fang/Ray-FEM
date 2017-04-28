function [ray, T, dTdx, dTdy] = eikonal_cgv(v0, G0, node0, node)
%% The analytical solution of eikonal equation with 
% constant gradient of velocity (CGV)
% 
%         v = v0 + 2*G0 \cdot (\x - \x0)
% 
% S02: 1x1
% g0: 1x2
% node0: 1x2
% node: Nx2

x = node(:,1);  y = node(:,2);
x0 = node0(:,1);  y0 = node0(:,2);  S0 = 1/v0;  
r = sqrt( (x-x0).^2 + (y-y0).^2 );       % |\x - \x0|
r2 = r.*r;
Gr = sqrt( G0(1)^2 + G0(2)^2 );          % |G0|
drdx = (x-x0)./r;  drdy = (y-y0)./r;     % dr/dx, dr/dy 

S = 1./( v0 + G0(1)*(x-x0) + G0(2)*(y-y0) );
dSdx = -S.*S*G0(1);  dSdy = -S.*S*G0(2);

Cz = 1/2*S0*Gr*Gr;
z = 1 + Cz*S.*r2;
dzdx = Cz*( dSdx.*r2 + 2*S.*r.*drdx);
dzdy = Cz*( dSdy.*r2 + 2*S.*r.*drdy);

T = 1/Gr*arccosh(z,0);
dTdx = 1/Gr*arccosh(z,1).*dzdx;
dTdy = 1/Gr*arccosh(z,1).*dzdy;
ray = atan2(dTdy, dTdx);  
ray = exp(1i*ray);
ray( r < 10*eps) = 0;


