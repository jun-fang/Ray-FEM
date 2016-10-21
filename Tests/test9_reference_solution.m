%% Reference solution for one point source problem: inside domain in
% heterogeneous medium

xs = 1/8;   ys = 1/10;             % point source location
speed = @(x) (3 - 2.5*exp( -((x(:,1)+1/8).^2 + (x(:,2)-0.1).^2)/0.8^2 )) ;    % wave-speed
cmin = 1/2;

a= 1/2;
plt = 0;                           % plot solution or not
fquadorder = 3; 
omega = 60*pi;

wpml = 2*(2*pi*cmin)/omega;        % width of PML
sigmaMax = 25/wpml;                % Maximun absorbtion

NPW = 40;
rh = 1/(NPW*round(omega/(2*pi*cmin)));
fprintf('omega/(2*pi) = %d,   1/h = %d,   NPW = %d \n',omega/(2*pi), 1/rh, NPW);
[node,elem] = squaremesh([-a,a,-a,a],rh);
[ru] = Standard_FEM_PML_PointSource(node,elem,omega,wpml,sigmaMax,xs,ys,speed,fquadorder,plt);
clear node elem;
save('test9_reference_solution.mat','rh','ru');
