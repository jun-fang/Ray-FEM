%% Test for accuracy of S-FEM with PML
clear;
fquadorder = 9;                  % numerical quadrature order
omega = 20*pi;              % high frequency
wl = 2*pi/omega;
wpml = 25/128;
sigmaMax = 25/wpml;
a = 1/2;

xs = 0; ys = 0;
% speed = @(x) ones(size(x(:,1)));
speed = @(x) 4/3*( 1-1/2*exp( -32*(x(:,1)-xs).^2 + (x(:,2)-ys).^2 ) );
% source = @(x) exp(-(4*omega/pi)^2*( (x(:,1)-xs).^2 + (x(:,2)-ys).^2));
nwl = 1/8;
sigma = nwl*wl;
source = @(x) 1/(2*pi*sigma^2)*exp(-( (x(:,1)-xs).^2 + (x(:,2)-ys).^2 )/2/sigma^2);

rpow = 11;
rh = 1/2^rpow;
[rnode,relem] = symmetric_squaremesh([-a,a,-a,a],rpow);
% rnode = [0,0; 1,0; 1,1; 0,1];
% relem = [2,3,1; 4,1,3];
% for k = 1:rpow
%     [rnode,relem] = uniformbisect(rnode,relem);
% end

fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('S-FEM:\n omega/(2*pi) = %.2d,   1/rh = %d   NPW = %.2d \n',omega/(2*pi), 1/rh, wl/rh);

tic;
rA = assemble_Helmholtz_matrix_SFEM(rnode,relem,omega,wpml,sigmaMax,speed,fquadorder);
rb = assemble_RHS(rnode,relem,source,fquadorder);
toc;
[~,~,isBdNode] = findboundary(relem);
freeNode = find(~isBdNode);
ru = zeros(size(rnode(:,1)));
tic;
ru(freeNode) = rA(freeNode,freeNode)\rb(freeNode);
toc;
[~,tru] = transform_node(rnode,ru);

mru = norm(ru,inf)

clear rA rb rnode relem freeNode isBdNode

n = 4;
rec_err = zeros(n,1);
rec_h = zeros(n,1);
for i = 1:n
    h = 1/2^(11-n-1+i);
    nh = 11-n-1+i;
    [node,elem] = symmetric_squaremesh([-a,a,-a,a],nh);
    fprintf(['\n' '-'*ones(1,80) '\n']);
    fprintf('S-FEM:\n omega/(2*pi) = %.2d,   1/rh = %d   NPW = %.2d \n',omega/(2*pi), 1/h, wl/h);
    
    tic;
    A = assemble_Helmholtz_matrix_SFEM(node,elem,omega,wpml,sigmaMax,speed,fquadorder);
    b = assemble_RHS(node,elem,source,fquadorder);
    toc;
    [~,~,isBdNode] = findboundary(elem);
    freeNode = find(~isBdNode);
    u = zeros(size(node(:,1)));
    tic;
    u(freeNode) = A(freeNode,freeNode)\b(freeNode);
    toc;
    [~,tu] = transform_node(node,u);
    norm(u,inf)
    
    rec_err(i) = getMAXerr(tu,h,tru,rh,wpml)/mru;
    rec_h(i) = h;
end

showrate(rec_h,rec_err);
fprintf('Mesh size 1/h  ');  fprintf('  &    %i   ',1./rec_h');  fprintf('\n');
fprintf('NPW            ');  fprintf('  &    %.1f  ',round(100*wl./rec_h)/100');  fprintf('\n');
fprintf('Max Error      ');  fprintf('  &  %.2e',rec_err*mru');  fprintf('\n');
fprintf('Rel Max Error  ');  fprintf('  &  %.2e',rec_err');  fprintf('\n');



