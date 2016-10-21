%% convergence rate of Ray-FEM
clear;
omegas = [4 8 16 32]'*2*pi;
no = length(omegas);
rfem_max_err = zeros(no,1);
rfem_l2_err = zeros(no,1);
sfem_max_err = zeros(no,1);
sfem_l2_err = zeros(no,1);


NPW = 16;
fquadorder = 9;                  % numerical quadrature order
wpml = 1/4;
sigmaMax = 25/wpml;
rh = 1/2048;
a = 1/2;
[rnode,relem] = squaremesh([-a,a,-a,a],rh);
    

xs = 0; ys = 0;
speed = @(x) ones(size(x(:,1)));
% source
% nwl = 1/16;  sigma = nwl*2*pi/omegas(end);
% sigma = rh/8;
% source = @(x) 1/(2*pi*sigma^2)*exp(-( (x(:,1)-xs).^2 + (x(:,2)-ys).^2 )/2/sigma^2);


starttime = cputime;
for ni = 1:no
    omega = omegas(ni);
    wl = 2*pi/omega;
    h = wl/NPW;
    [node,elem] = squaremesh([-a,a,-a,a],h);
    
    sigma = 1/2*h;
    dr = 4*sigma;
    source = @(x) 1/(2*pi*sigma^2)*exp(-( (x(:,1)-xs).^2 + (x(:,2)-ys).^2 )/2/sigma^2);
    fprintf(['\n' '-'*ones(1,80) '\n']);
    fprintf('Gaussian point source parameter sigma = 1/%.0f wavelength\n', wl/sigma);

    
    fprintf(['\n' '-'*ones(1,80) '\n']);
    fprintf('Ray-FEM:\n omega/(2*pi) = %.2d,   1/h = %d   NPW = %.2d \n',omega/(2*pi), 1/h, wl/h);

    
    %% Exact Ray-FEM
    tic;
    ray_ang = ex_ray(node,xs,ys);
    d = sqrt((node(:,1)-xs).^2 + (node(:,2)-ys).^2);
    ray = exp(1i*ray_ang).*(d>10*eps);
    A = assemble_Helmholtz_matrix_with_ray_1(node,elem,omega,wpml,sigmaMax,speed,ray,fquadorder);
    b = assemble_RHS_with_ray_1(node,elem,omega,source,speed,ray,fquadorder);
    fb = b;
    idx = find(d<dr);
    for i = 1:length(idx)
        ii = idx(i);
        xi = node(ii,1); yi = node(ii,2);
        fb(ii) = RHS_integral_with_ray(xi,yi,h,rh/16,omega,ray(ii),source,fquadorder);
    end
    norm(b-fb,inf)
    toc;
    [~,~,isBdNode] = findboundary(elem);
    freeNode = find(~isBdNode);
    v = zeros(size(node(:,1)));
    tic;
    v(freeNode) = A(freeNode,freeNode)\fb(freeNode);
    toc;
    
    % reconstruct the solution
    grad = ray(:);
    grad = [real(grad),imag(grad)];
    temp = grad(:,1).*node(:,1) + grad(:,2).*node(:,2);
    
    c = speed(node);    % medium speed
    k = omega./c;           % wavenumber
    u = v.*exp(1i*k(:).*temp);
    u = sum(u,2);
    
    %% S-FEM with same mesh size
    sA = assemble_Helmholtz_matrix_SFEM(node,elem,omega,wpml,sigmaMax,speed,fquadorder);
    sb = assemble_RHS(node,elem,source,fquadorder);
    fsb = sb;
    for i = 1:length(idx)
        ii = idx(i);
        xi = node(ii,1); yi = node(ii,2);
        fsb(ii) = RHS_integral(xi,yi,h,rh/16,source,fquadorder);
    end
    norm(sb-fsb,inf)
    su = zeros(size(node(:,1)));
    su(freeNode) = sA(freeNode,freeNode)\fsb(freeNode);
    
    
    %% Reference solution S-FEM
    tic;
    rA = assemble_Helmholtz_matrix_SFEM(rnode,relem,omega,wpml,sigmaMax,speed,fquadorder);
    rb = assemble_RHS(rnode,relem,source,fquadorder);
    toc;
    frb = rb;
    idx = find(sqrt((rnode(:,1)-xs).^2 + (rnode(:,2)-ys).^2)<dr);
    for i = 1:length(idx)
        ii = idx(i);
        xi = rnode(ii,1); yi = rnode(ii,2);
        fb(ii) = RHS_integral(xi,yi,rh,rh/16,source,fquadorder);
    end
    norm(rb-frb,inf)
    
    [~,~,isBdNode] = findboundary(relem);
    freeNode = find(~isBdNode);
    ru = zeros(size(rnode(:,1)));
    tic;
    ru(freeNode) = rA(freeNode,freeNode)\frb(freeNode);
    toc;
    
    rfem_max_err(ni) = getMAXerr(u,h,ru,rh,wpml)
    rfem_l2_err(ni) = getRelL2err(u,h,ru,rh,wpml)
    sfem_max_err(ni) = getMAXerr(su,h,ru,rh,wpml)
    sfem_l2_err(ni) = getRelL2err(su,h,ru,rh,wpml)
    
end

figure(10);    % convergence rate
subplot(2,2,1);
showrate(omegas,rfem_max_err);
subplot(2,2,2);
showrate(omegas,rfem_l2_err);
subplot(2,2,3);
showrate(omegas,sfem_max_err);
subplot(2,2,4);
showrate(omegas,sfem_l2_err);

endtime = cputime;

endtime - starttime


fprintf('R-FEM Max Error  ');  fprintf('  &  %.2e',rfem_max_err');  fprintf('\n');
fprintf('S-FEM Max Error  ');  fprintf('  &  %.2e',sfem_max_err');  fprintf('\n');
fprintf('R-FEM L2 Error   ');  fprintf('  &  %.2e',rfem_l2_err');  fprintf('\n');
fprintf('S-FEM L2 Error   ');  fprintf('  &  %.2e',sfem_l2_err');  fprintf('\n');
