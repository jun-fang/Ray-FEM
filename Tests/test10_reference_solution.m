%% Reference solution for one point source problem: inside domain in
% heterogeneous medium: caustics example

xs = -0.2;   ys = -0.3;             % point source location
speed = @(x) ( 1+ 0.5*sin(2*pi*x(:,1)));
cmin = 1/2;

a= 1;
plt = 0;                           % plot solution or not
fquadorder = 3; 
omega = 60*pi;

wpml = 10*(2*pi*cmin)/omega;        % width of PML
sigmaMax = 50/wpml;                % Maximun absorbtion

NPW = 15;
rh = 1/(NPW*round(omega/(2*pi*cmin)));
fprintf('omega/(2*pi) = %d,   1/h = %d,   NPW = %d \n',omega/(2*pi), 1/rh, NPW);
[rnode,relem] = squaremesh([-a,a,-a,a],rh);
[ru] = Standard_FEM_PML_PointSource(rnode,relem,omega,wpml,sigmaMax,xs,ys,speed,fquadorder,plt);
% clear node elem;
save('test10_reference_solution_20_1.mat','rh','ru');
% test10_caustics;



