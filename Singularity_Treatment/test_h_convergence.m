global omega a;
xs = 0;   ys = 0;
omega = 10*pi;
speed = @(p) ones(size(p,1),1);    % medium speed


% exu = @(p) ( exp(p(:,1))+exp(p(:,2)) )...
%     .*exp( 1i*omega*sqrt(2)/2*( p(:,1)+p(:,2) ) );
% source = @(p) -(1+1i*sqrt(2)*omega).*exu(p);

fquadorder = 3;                    % numerical quadrature order
pde = Helmholtz_data9;
solver = 'DIR';
plt = 0;

omega = 10*pi;
epsilon = 0.131;
wl = 2*pi/omega;

NPW = 6;
h = wl/NPW;
h = 1/round(1/h);

wpml = 0.2;%3*round(wavelength/h)*h;                 % width of PML
sigmaMax = 25/wpml;                % Maximun absorbtion

a = 1;

nt = 3;

hs = zeros(1,nt);
error = zeros(4,nt);
for ii = 1:nt
    ii
    tic;
    h = h/2;
    hs(ii) = h;
    [node,elem] = squaremesh([-a,a,-a,a],h);
    N = size(node,1);
%     [u,A,b,rel_L2_err] = Standard_FEM_IBC(node,elem,omega,pde,fquadorder,solver,plt);
%     error(4,ii) = rel_L2_err;

    [u,A,b] = Standard_FEM_PML(node,elem,omega,wpml,sigmaMax,pde,fquadorder,solver,plt);
    
    
%     A = assemble_Helmholtz_matrix_SFEM(node,elem,omega,wpml,sigmaMax,speed,fquadorder);
%     b = assemble_RHS_PML(node,elem,omega,wpml,sigmaMax,source,fquadorder);
    
%     %% Boundary conditions
%     [bdNode,~,isBdNode] = findboundary(elem);
%     freeNode = find(~isBdNode);
%     u = zeros(N,1);
%     u(freeNode) = A(freeNode,freeNode)\b(freeNode);
%     
    v = exu(node);
    
    du = u-v;
    x = node(:,1); y = node(:,2);
    du(x>=max(x)-wpml)=0; du(x<= min(x)+wpml) = 0;
    du(y>=max(y)-wpml)=0; du(y<= min(y)+wpml) = 0;
    
    max_err = norm(du,inf);
    rel_max_err = norm(du,inf)/norm(v,inf);
    l2_err = norm(du)*h;
    rel_l2_err = norm(du)/norm(v);
    
    error(:,ii) = [max_err,rel_max_err,l2_err,rel_l2_err]';
    toc;
end

error

subplot(2,2,1);
FJ_showrate(hs,error(1,:));
subplot(2,2,2);
FJ_showrate(hs,error(2,:));
subplot(2,2,3);
FJ_showrate(hs,error(3,:));
subplot(2,2,4);
FJ_showrate(hs,error(4,:));




