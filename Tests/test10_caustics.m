xs = -0.2;   ys = -0.3;             % point source location
speed = @(x) ( 1+ 0.4*sin(4*pi*(x(:,1) + 0)).*sin(2*pi*x(:,2)) );
% speed = @(x) (1 - 0.5*exp(-(x(:,1)-1).^2));
speed = @(x) ( 1+ 0.5*sin(2*pi*(x(:,1) + 0.0)));


a= 1/2;
plt = 0;                           % plot solution or not
fquadorder = 3; 
omega = 60*pi;
wpml = 2*pi/omega;                 % width of PML
sigmaMax = 25/wpml;                % Maximun absorbtion
cmin = 1/2;

NPW = 10;
h = 1/(NPW*round(omega/(2*pi*cmin)));
1/h
[node,elem] = squaremesh([-a,a,-a,a],h);
% [node,elem] = squaremesh([0,1,0,2],h);
[u] = Standard_FEM_PML_PointSource(node,elem,omega,wpml,sigmaMax,xs,ys,speed,fquadorder,plt);

% clear node elem;

[X,Y] = meshgrid(-a:h:a,-a:h:a);
% [X,Y] = meshgrid(0:h:1,0:h:2);
[m,n] = size(X);
uh = reshape(u,m,n);

figure(1)
contour(X,Y,real(uh));
title('Level set of reference solution u_{ref}');
axis equal;

figure(2);
FJ_showresult(node,elem,real(u));
axis equal;

% save('test9_reference_solution.mat','u','h');
