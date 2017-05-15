function [theta,r,uu] = mapto_polar(node,elem,omega,speed,v,ray,xs,ys,r1,r2,theta1,theta2,u,method)

[theta,r] = meshgrid(theta1:1/10000:theta2, r1: 1/5000:r2);
[m,n] = size(r);
xx = r.*cos(theta) + xs;
yy = r.*sin(theta) + ys;
xynode = [xx(:),yy(:)];

if nargin == 12
    uu = RayFEM_solution(node,elem,omega,speed,v,ray,xynode);
else
    xmin = min(node(:,1)); xmax = max(node(:,1));
    ymin = min(node(:,2)); ymax = max(node(:,2));
    h = (xmax-xmin)/(round(sqrt(length(node(:,1))))-1);
    [X,Y] = meshgrid(xmin:h:xmax,ymin:h:ymax);
    [Xm,Xn] = size(X);
    V = reshape(u,Xm,Xn);
    uu = interp2(X,Y,V,xx(:),yy(:),method);
end

uu = reshape(uu,m,n);



mesh(theta,r,real(uu));
xlabel('theta');
ylabel('r');
% axis equal;
axis tight;
az = 0;
el = 90;
view(az, el);

% colorbar;

