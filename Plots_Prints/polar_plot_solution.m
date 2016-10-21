function [theta,rn,uu] = polar_plot_solution(omega,node,xs,ys,r1,r2,theta1,theta2,u,method,opt)

wl = 2*pi/omega;
[theta,r] = meshgrid(theta1:1/10000:theta2, r1: 1/25000:r2);
[m,n] = size(r);
xx = r.*cos(theta) + xs;
yy = r.*sin(theta) + ys;
% xynode = [xx(:),yy(:)];
rn = r/wl;

xmin = min(node(:,1)); xmax = max(node(:,1));
ymin = min(node(:,2)); ymax = max(node(:,2));
h = (xmax-xmin)/(round(sqrt(length(node(:,1))))-1);
[X,Y] = meshgrid(xmin:h:xmax,ymin:h:ymax);
[Xm,Xn] = size(X);
V = reshape(u,Xm,Xn);
uu = interp2(X,Y,V,xx(:),yy(:),method);
uu = reshape(uu,m,n);


if opt == 1
    mesh(theta/pi,rn,real(uu));
    xlabel('\theta/\pi');
    ylabel('r/\lambda');
else 
    mesh(theta/pi,rn*wl,real(uu));
    xlabel('\theta/\pi');
    ylabel('r');
end
% axis equal;
axis tight;
az = 0;
el = 90;
view(az, el);

% colorbar;

