

xs = 0; ys = 0;
omega = 20*pi;
wl = 2*pi/omega;
sigma = wl/8;
speed = @(p) ones(size(p(:,1)));
source = @(p) 1/(2*pi*sigma*sigma)...
    *exp( -( (p(:,1)-xs).^2 + (p(:,2)-ys).^2 )/(2*sigma*sigma) );

wpml = 25/128;
sigmaMax = 25/wpml;
fquadorder = 3;
a = 1/2;
load('test_PML.mat');


%% get reference solution
% tic;
% hr = 1/2048;
% [node,elem] = squaremesh([-a,a,-a,a],hr);
% A = assemble_Helmholtz_matrix_SFEM(node,elem,omega,wpml,sigmaMax,speed,fquadorder);
% b = assemble_RHS_PML(node,elem,omega,wpml,sigmaMax,source,fquadorder);
% 
% %% Boundary conditions
% [bdNode,~,isBdNode] = findboundary(elem);
% freeNode = find(~isBdNode);
% N = size(node,1);  n = round(sqrt(N));
% u = zeros(N,1);
% u(freeNode) = A(freeNode,freeNode)\b(freeNode);
% ur = reshape(u,n,n);
% url2 = norm(u)*hr;
% toc;
% 
% save('test_PML.mat','hr','ur','url2');
% 
% clear node elem u;


%% test h convergence 
nt = 4;
hs = zeros(nt,1);
errs_com = zeros(nt,1);  % error on computational domain
errs_phy = zeros(nt,1);  % error on physical domain
h = 1/64;
for ii = 1:nt
    ii
    tic;
    h = h/2;
    hs(ii) = h;
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
    idx = h/hr*(idx-1)+1;
    du = u - ur(idx,idx);  
    du = du(:);
    errs_com(ii) = norm(du)*h;
    
    x = node(:,1); y = node(:,1);
    idx = find( (x<=max(x)-wpml).*(x>= min(x)+wpml)...
        .*(y>= min(y)+wpml).*(y<=max(y)+wpml) ); 
    du_phy = du(idx);
    errs_phy(ii) = norm(du_phy)*h;
    toc;
end

figure(10);
subplot(1,2,1);
showrate(hs,errs_com);

subplot(1,2,2);
showrate(hs,errs_phy);

    
