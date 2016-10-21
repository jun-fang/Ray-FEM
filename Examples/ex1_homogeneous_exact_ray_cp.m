%% Test for Ray-FEM in homogeneous medium 
% Gausssian point source with exact ray information
% 1. Gaussian parameter sigma should scale with frequency to approximate it
% as a point source
% 2. High order numerical quadrature rule or very fine mesh to compute the
% integral near the source location 
%
% 

clear;
fquadorder = 9;                  % numerical quadrature order
omega = 20*pi;              % high frequency
wl = 2*pi/omega;
wpml = 25/128;
sigmaMax = 25/wpml;
a = 1/2;
 
r1 = 0; r2 = 0;
speed = @(x) ones(size(x(:,1)));
% speed = @(x) 4/3*( 1-1/2*exp( -32*(x(:,1)-r1).^2 + (x(:,2)-r2).^2 ) );
nwl = 1/16;
sigma = nwl*wl;
source = @(x) 1/(2*pi*sigma^2)*exp(-( (x(:,1)-r1).^2 + (x(:,2)-r2).^2 )/2/sigma^2);

rpow = 11;
rh = 1/2^rpow;
[rnode,relem] = symmetric_squaremesh([-a,a,-a,a],rpow);
% rnode = [0,0; 1,0; 1,1; 0,1];
% relem = [2,3,1; 4,1,3];
% for k = 1:rpow
%     [rnode,relem] = uniformbisect(rnode,relem);
% end

fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('Gaussian point source parameter sigma = %f wavelength\n', sigma/wl);

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
rec_ray_max_err = zeros(n,1); rec_ray_l2_err = zeros(n,1);
rec_s_max_err = zeros(n,1);   rec_s_l2_err = zeros(n,1);
rec_h = zeros(n,1);
rec_ray_u = cell(n,1);  
rec_s_u = cell(n,1);

for i = 1:n
    h = 1/2^(rpow-n-1+i);
    hp = rpow-n-1+i;
    [node,elem] = symmetric_squaremesh([-a,a,-a,a],hp);
    fprintf(['\n' '-'*ones(1,80) '\n']);
    fprintf('Ray-FEM:\n omega/(2*pi) = %.2d,   1/rh = %d   NPW = %.2d \n',omega/(2*pi), 1/h, wl/h);
        
    
    %% Exact Ray-FEM
    tic;
    ray_ang = ex_ray(node,r1,r2);
    d = sqrt((node(:,1)-r1).^2 + (node(:,2)-r2).^2);
    ray = exp(1i*ray_ang).*(d>10*eps);
    A = assemble_Helmholtz_matrix_with_ray_1(node,elem,omega,wpml,sigmaMax,speed,ray,fquadorder);
    b = assemble_RHS_with_ray_1(node,elem,omega,source,speed,ray,fquadorder); 
    % high order or accurate numerical quadrature to compute the right hand
    % side
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
    [~,tu] = transform_node(node,u);
    rec_ray_u{i} = reshape(tu,round(sqrt(length(tu))),round(sqrt(length(tu))));

    fprintf('Ray-FEM max norm: %f \n', norm(u,inf));
    
    
    
    %% S-FEM
    sA = assemble_Helmholtz_matrix_SFEM(node,elem,omega,wpml,sigmaMax,speed,fquadorder);
    sb = assemble_RHS(node,elem,source,fquadorder);
    su = zeros(size(node(:,1)));
    su(freeNode) = sA(freeNode,freeNode)\sb(freeNode);
    [~,tsu] = transform_node(node,su);
    rec_s_u{i} = reshape(tsu,round(sqrt(length(tsu))),round(sqrt(length(tsu))));

    fprintf('S-FEM max norm: %f \n', norm(su,inf));
    
    
    rec_ray_max_err(i) = getMAXerr(tu,h,tru,rh,wpml);
    rec_s_max_err(i) = getMAXerr(tsu,h,tru,rh,wpml);
    rec_ray_l2_err(i) = getRelL2err(tu,h,tru,rh,wpml);
    rec_s_l2_err(i) = getRelL2err(tsu,h,tru,rh,wpml);
    rec_h(i) = h;
end

%% print out results
fprintf('Mesh size 1/h   ');  fprintf('  &    %i   ',1./rec_h');  fprintf('\n');
fprintf('NPW             ');  fprintf('  &    %.1f  ',round(100*wl./rec_h)/100');  fprintf('\n');
fprintf('R-FEM Max Error  ');  fprintf('  &  %.2e',rec_ray_max_err');  fprintf('\n');
fprintf('S-FEM Max Error  ');  fprintf('  &  %.2e',rec_s_max_err');  fprintf('\n');
fprintf('R-FEM L2 Error   ');  fprintf('  &  %.2e',rec_ray_l2_err');  fprintf('\n');
fprintf('S-FEM L2 Error   ');  fprintf('  &  %.2e',rec_s_l2_err');  fprintf('\n');


%% plot
figure(10);    % convergence rate
subplot(2,2,1);
showrate(rec_h,rec_ray_max_err);
subplot(2,2,2);
showrate(rec_h,rec_ray_l2_err);
subplot(2,2,3);
showrate(rec_h,rec_s_max_err);
subplot(2,2,4);
showrate(rec_h,rec_s_l2_err);



figure(20);   % wave-field
yy = 0.29;
yy = rec_h(1)*ceil((yy+a)/rec_h(1));

% reference solution
rn = round(sqrt(length(tru)));
ref_u = reshape(tru,rn,rn);
ryn = round(yy/rh) + 1;
rxx = -a:rh:a;
ruh = ref_u(ryn,:);

plot(rxx,real(ruh),'k');
hold on;

% solution 1
yn = round(yy/rec_h(1)) + 1;
xx = -a:rec_h(1):a;
uh = rec_ray_u{1};
uh = uh(yn,:);
plot(xx,real(uh),'gx-');
hold on;
uh = rec_s_u{1};
uh = uh(yn,:);
plot(xx,real(uh),'g*:');
hold on;

% solution 2
yn = round(yy/rec_h(2)) + 1;
xx = -a:rec_h(2):a;
uh = rec_ray_u{2};
uh = uh(yn,:);
plot(xx,real(uh),'mo-');
hold on;
uh = rec_s_u{2};
uh = uh(yn,:);
plot(xx,real(uh),'md:');
hold on;

% solution 3
yn = round(yy/rec_h(3)) + 1;
xx = -a:rec_h(3):a;
uh = rec_ray_u{3};
uh = uh(yn,:);
plot(xx,real(uh),'bs-');
hold on;
uh = rec_s_u{3};
uh = uh(yn,:);
plot(xx,real(uh),'b^:');
hold on;

% solution 4
yn = round(yy/rec_h(4)) + 1;
xx = -a:rec_h(4):a;
uh = rec_ray_u{4};
uh = uh(yn,:);
plot(xx,real(uh),'r+-');
hold on;
uh = rec_s_u{4};
uh = uh(yn,:);
plot(xx,real(uh),'r<:');
hold on;

xlabel('x');
ylabel('Real part wavefield');
legend('Reference solution','Ray-FEM h_1','S-FEM h_1',...
    'Ray-FEM h_2','S-FEM h_2','Ray-FEM h_3','S-FEM h_3',...
    'Ray-FEM h_4','S-FEM h_4','LOCATION','Best');
title(['Wavefield at y = ' num2str(yy-a)],'FontSize', 14)

