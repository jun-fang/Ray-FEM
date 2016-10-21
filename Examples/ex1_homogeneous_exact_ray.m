%% Test for accuracy of S-FEM with PML
clear;
fquadorder = 9;                  % numerical quadrature order
omega = 20*pi;              % high frequency
wl = 2*pi/omega;
wpml = 1/4; %25/128;
sigmaMax = 25/wpml;
a = 1/2;

r1 = 0; r2 = 0;
speed = @(x) ones(size(x(:,1)));
% speed = @(x) 4/3*( 1-1/2*exp( -32*(x(:,1)-r1).^2 + (x(:,2)-r2).^2 ) );
source = @(x) exp(-(omega/pi/4)^2*( (x(:,1)-r1).^2 + (x(:,2)-r2).^2));

rpow = 12;
rh = 1/2^rpow;
[rnode,relem] = squaremesh([-a,a,-a,a],rh);
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

norm(ru,inf)


% ray_ang = ex_ray(rnode,r1,r2);
% d = sqrt((rnode(:,1)-r1).^2 + (rnode(:,2)-r2).^2);
% ray = exp(1i*ray_ang).*(d>10*eps);
% rA = assemble_Helmholtz_matrix_with_ray_1(rnode,relem,omega,wpml,sigmaMax,speed,ray,fquadorder);
% rb = assemble_RHS_with_ray_1(rnode,relem,omega,source,speed,ray,fquadorder);
% toc;
% 
% [~,~,isBdNode] = findboundary(relem);
% freeNode = find(~isBdNode);
% rv = zeros(size(rnode(:,1)));
% tic;
% rv(freeNode) = rA(freeNode,freeNode)\rb(freeNode);
% toc;
% 
% %% reconstruct the solution
% grad = ray(:);
% grad = [real(grad),imag(grad)];
% temp = grad(:,1).*rnode(:,1) + grad(:,2).*rnode(:,2);
% 
% c = speed(rnode);    % medium speed
% k = omega./c;           % wavenumber
% 
% ru = rv.*exp(1i*k(:).*temp);


clear rA rb rnode relem freeNode isBdNode

rpow = 10;
n = 4;
rec_err = zeros(n,1);
rec_h = zeros(n,1);
for i = 1:n
    h = 1/2^(rpow-n-1+i);
    [node,elem] = squaremesh([-a,a,-a,a],h);
    fprintf(['\n' '-'*ones(1,80) '\n']);
    fprintf('Ray-FEM:\n omega/(2*pi) = %.2d,   1/rh = %d   NPW = %.2d \n',omega/(2*pi), 1/h, wl/h);
    
    
    %% Exact Ray-FEM
    tic;
    ray_ang = ex_ray(node,r1,r2);
    d = sqrt((node(:,1)-r1).^2 + (node(:,2)-r2).^2);
    ray = exp(1i*ray_ang).*(d>10*eps);
    A = assemble_Helmholtz_matrix_with_ray_1(node,elem,omega,wpml,sigmaMax,speed,ray,fquadorder);
    b = assemble_RHS_with_ray_1(node,elem,omega,source,speed,ray,fquadorder);
    toc;
    [~,~,isBdNode] = findboundary(elem);
    freeNode = find(~isBdNode);
    v = zeros(size(node(:,1)));
    tic;
    v(freeNode) = A(freeNode,freeNode)\b(freeNode);
    toc;
    
    %% reconstruct the solution
    grad = ray(:);
    grad = [real(grad),imag(grad)];
    temp = grad(:,1).*node(:,1) + grad(:,2).*node(:,2);
    
    c = speed(node);    % medium speed
    k = omega./c;           % wavenumber
    
    u = v.*exp(1i*k(:).*temp);
    u = sum(u,2);
    norm(u,inf)
    
    rec_err(i) = getMAXerr(u,h,ru,rh);
    rec_h(i) = h;
end

showrate(rec_h,rec_err);
fprintf('Mesh size 1/h  ');  fprintf('  &    %i   ',1./rec_h');  fprintf('\n');
fprintf('NPW            ');  fprintf('  &    %.1f  ',round(100*wl./rec_h)/100');  fprintf('\n');
fprintf('Max Error      ');  fprintf('  &  %.2e',rec_err');  fprintf('\n');



