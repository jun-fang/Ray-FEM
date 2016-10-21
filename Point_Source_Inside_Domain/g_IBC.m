function g = g_IBC(omega,p,xc,yc,x1,x2,y1,y2)
%% ???? normal vector?????

k = omega;
xx = p(:,1); yy = p(:,2);
u = sqrt(k)*besselh(0,1,k*sqrt((xx-xc).^2 + (yy-yc).^2));

x = p(:,1)-xc; y = p(:,2)-yc;
r = sqrt(x.^2 + y.^2);

Dux = -sqrt(k)*besselh(1,1,k*r).*k.*x./r;
Duy = -sqrt(k)*besselh(1,1,k*r).*k.*y./r;

uN = Dux;

east_idx = find(abs(xx-x2)<10*eps);
uN(east_idx) = -Dux(east_idx);
west_idx = find(abs(xx-x1)<10*eps);
uN(west_idx) = Dux(west_idx);

north_idx = find(abs(yy-y2)<10*eps);
uN(north_idx) = -Duy(north_idx);
south_idx = find(abs(yy-y1)<10*eps);
uN(south_idx) = Duy(south_idx);

g = uN + 1i*k*u;