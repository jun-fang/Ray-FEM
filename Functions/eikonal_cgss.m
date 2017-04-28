function [ray, T, dT1dx, dT1dy, dT2dx, dT2dy] = eikonal_cgss(S02, g0, node0, node)
%% The analytical solution of eikonal equation with 
% constant gradient of slowness squared(CGSS)
% 
%         S^2 = |\nabla T|^2 = S0^2 + 2*g0 \cdot (\x - \x0)
% 
% two branch of solutions: T(\x) = S_bar^2 * sigma - 1/6 * |g0|^2 * sigma^3
% 
% S02: 1x1
% g0: 1x2
% node0: 1x2
% node: Nx2

x = node(:,1);  y = node(:,2);
x0 = node0(:,1);  y0 = node0(:,2);
r = sqrt( (x-x0).^2 + (y-y0).^2 );       % |\x - \x0|
drdx = (x-x0)./r;  drdy = (y-y0)./r;   % dr/dx, dr/dy 

S_bar = sqrt( S02 + g0(1)*(x-x0) + g0(2)*(y-y0) );
dSdx = g0(1) ./ (2*S_bar);  dSdy = g0(2) ./ (2*S_bar);

gr = sqrt( g0(1)^2 + g0(2)^2 );
D = sqrt( S_bar.^4 - gr^2*r.^2 );
dDdx = ( 2*S_bar.^3.*dSdx - gr^2*r.*drdx )./D;
dDdy = ( 2*S_bar.^3.*dSdy - gr^2*r.*drdy )./D;


%% branch one 
s1 = sqrt( 2*(S_bar.^2 + D) );
ds1dx = ( 2*S_bar.*dSdx + dDdx )./s1;
ds1dy = ( 2*S_bar.*dSdy + dDdy )./s1;

T1 = ( S_bar.^2.*s1 - s1.^3/6 )/gr;
dT1dx = ( 2*S_bar.*s1.*dSdx + S_bar.^2.*ds1dx - s1.^2.*ds1dx/2 )/gr;
dT1dy = ( 2*S_bar.*s1.*dSdy + S_bar.^2.*ds1dy - s1.^2.*ds1dy/2 )/gr;


%% branch two 
s2 = 2*r./s1; % sqrt( 2*(S_bar.^2 - D) );
ds2dx = 2*( s1.*drdx - r.*ds1dx )./(s1.*s1);
ds2dy = 2*( s1.*drdy - r.*ds1dy )./(s1.*s1);

T2 = s2.*( S_bar.^2 - 1/6*(gr*s2).^2);
dT2dx = 2*S_bar.*s2.*dSdx + ( S_bar.^2 - 1/2*(gr*s2).^2 ).*ds2dx;
dT2dy = 2*S_bar.*s2.*dSdy + ( S_bar.^2 - 1/2*(gr*s2).^2 ).*ds2dy;


% T2 = ( S_bar.^2.*s2 - s2.^3/6 )/gr;
% dT2dx = ( 2*S_bar.*s2.*dSdx + S_bar.^2.*ds2dx - s2.^2.*ds2dx/2 )/gr;
% dT2dy = ( 2*S_bar.*s2.*dSdy + S_bar.^2.*ds2dy - s2.^2.*ds2dy/2 )/gr;


%% two braches:
T = [T1, T2];
ray1 = atan2(dT1dy, dT1dx);  
ray1( r < 10*eps) = 0;
ray1 = exp(1i*ray1);
ray2 = atan2(dT2dy, dT2dx);  ray2 = exp(1i*ray2);
ray1( r < 10*eps) = 0;
ray2( r < 10*eps) = 0;
ray = [ray1, ray2]; 



