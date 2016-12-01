xs = 0; y0 = 0;
omega = 20*pi;
wl = 2*pi/omega;
sigma = wl/8;
source = @(p) 1/(2*pi*sigma*sigma)...
    *exp( -( (p(:,1)-xs).^2 + (p(:,2)-ys).^2 )/(2*sigma*sigma) );

wpml = 25/128;
sigmaMax = 25/wpml;
fquadorder = 3;
a = 1/2;

tic
hr = 1/2048;
[node,elem] = squaremesh([-a,a,-a,a],hr);
A = assemble_Helmholtz_matrix_SFEM(node,elem,omega,wpml,sigmaMax,speed,fquadorder);
b = assemble_RHS_PML(node,elem,omega,wpml,sigmaMax,source,fquadorder);

%% Boundary conditions
[bdNode,~,isBdNode] = findboundary(elem);
freeNode = find(~isBdNode);
N = size(node,1);  n = round(sqrt(N));
u = zeros(N,1);
u(freeNode) = A(freeNode,freeNode)\b(freeNode);
ur = reshape(u,n,n);
toc;

clear node elem ur;

nt = 4;
hs = zeros(nt,1);
h = 1/64;
for ii = 1:nt
    ii
    h = h/2;
    [node,elem] = squaremesh([-a,a,-a,a],h);
    A = assemble_Helmholtz_matrix_SFEM(node,elem,omega,wpml,sigmaMax,speed,fquadorder);
    b = assemble_RHS_PML(node,elem,omega,wpml,sigmaMax,source,fquadorder);
    
    %% Boundary conditions
    [bdNode,~,isBdNode] = findboundary(elem);
    freeNode = find(~isBdNode);
    N = size(node,1);  n = round(sqrt(N));
    u = zeros(N,1);
    u(freeNode) = A(freeNode,freeNode)\b(freeNode);
    u = reshape(u,n,n);
    idx = 1:n;
    idx = h/rh*(idx-1)+1;
    
    
