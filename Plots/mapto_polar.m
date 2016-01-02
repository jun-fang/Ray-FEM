function [theta,r,uu] = mapto_polar(node,elem,omega,speed,v,ray,xs,ys,r1,r2,theta1,theta2)

[theta,r] = meshgrid(theta1:1/10000:theta2, r1: 1/5000:r2);
[m,n] = size(r);
xx = r.*cos(theta) + xs;
yy = r.*sin(theta) + ys;
xynode = [xx(:),yy(:)];
uu = ray_solution(node,elem,omega,speed,v,ray,xynode);
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

