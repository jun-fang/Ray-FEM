%% A test example for Ray-FEM: assembling with ray information

rec_omega = [ 25 50 100 200 400]*pi;
high_omega = rec_omega;
n = 5;
max_err = zeros(1,n);
rel_max_err = zeros(1,n);
l2_err = zeros(1,n);
rel_l2_err = zeros(1,n);
for i = 1:n
    i
    tic;
    
    omega = rec_omega(i);
    xc = 0; yc = 0;
    sigma = 1/100;
    
    %% Load source data
    u_exact = @(x) besselh(0,1,omega*sqrt((x(:,1)-xc).^2 + (x(:,2)-yc).^2));
    %u_pl = @(x) (2/pi./(sqrt((x(:,1)-xc).^2 + (x(:,2)-yc).^2))).^0.5.*exp(1i*(omega*sqrt((x(:,1)-xc).^2 + (x(:,2)-yc).^2) - pi/4));
%     source = @(x) -4*sqrt(omega)*1i*1/(2*pi*sigma^2)*exp( -( (x(:,1)-xc).^2 + (x(:,2)-yc).^2  )/(2*sigma^2) );
    source = @(x) 0*ones(size(x(:,1)));
    speed = @(x) ones(size(x(:,1)));
    
    %% Set up
    plt = 0;                   % show solution or not
    fquadorder = 3;            % numerical quadrature order
    NPW = 8;
    h = 1/round((NPW*omega)/(2*pi));
    
    wavelength = 2*pi/omega;
    wpml = 4*round(wavelength/h)*h;               % width of PML
    sigmaMax = 50/wpml;                % Maximun absorbtion
    
    low_omega = 2*sqrt(omega);
    bd = 1/omega^(1/4);
    bd = round(bd/h)*h;
    bd = 0.3;
    wpml = 0.1;  bd = 0.2;
    
    a = 1/2+wpml;
    [node,elem] = squaremesh_annulus([-a,a,-a,a],[-bd,bd,-bd,bd],h);
    
    
    %% Exact ray information
    xx = node(:,1)-xc; yy = node(:,2)-yc;
    ray = atan2(yy,xx);
    ray = exp(1i*ray);
    ray = ray.*(1 - (abs(xx)<eps).*(abs(yy)<eps));
    
    N = size(node,1);       % number of grid points
    Nray = size(ray,2);     % number of rays crossing at each grid node
    Ndof = N*Nray;
    
    
    %% Assembling the linear system    
    A = assemble_Helmholtz_matrix_with_ray(node,elem,omega,wpml,sigmaMax,speed,ray,fquadorder);
    b = assemble_RHS_with_ray(node,elem,omega,source,speed,ray,fquadorder);
    % [A,b] = modify_linear_system_with_boundary_conditons(A,b,xc,yc,node,elem,omega,speed,ray,fquadorder,'Neu');
    
    
    %% Boundaries
    [~,bdEdge,isBdNode] = findboundary(elem);
    xmax = max(node(:,1)); xmin = min(node(:,1));
    ymax = max(node(:,2)); ymin = min(node(:,2));
    
    bx = node(:,1); by = node(:,2);
    isBdNode1 = ( bx>xmin+10*eps ).*( bx<xmax-10*eps )...
        .*( by>ymin+10*eps ).*( by<ymax-10*eps );
    isBdNode1 = 1 - isBdNode1;
    x = node(:,1); y = node(:,2);
    isBdNode2 = (x>=-bd-10*eps).*(x<=bd+10*eps).*(y>=-bd-10*eps).*(y<=bd+10*eps);
    
    rep_isBdNode = repmat(isBdNode,1,Nray);
    isBdNode = rep_isBdNode(:);
    freeNode = find(~isBdNode);
    
    idx2 = isBdNode2>0;
    BdNode2 = node(idx2,:);
    dx = real(ray(idx2)).*BdNode2(:,1) + imag(ray(idx2)).*BdNode2(:,2);
    
    
    %% Solve Av=b and reconstruct the solution
    v = zeros(Ndof,1);
    v(idx2) = u_exact(BdNode2).*exp(-1i*omega*dx);
    b = b - A*v;
    v(freeNode) = A(freeNode,freeNode)\b(freeNode);
    
    grad = ray(:);
    grad = [real(grad),imag(grad)];
    repnode = repmat(node,Nray,1);
    temp = grad(:,1).*repnode(:,1) + grad(:,2).*repnode(:,2);
    
    c = speed(node);    % medium speed
    k = omega./c;           % wavenumber
    kk = repmat(k,1,Nray);
    
    u = v.*exp(1i*kk(:).*temp);
    u = reshape(u,N,Nray);
    u = sum(u,2);
    
    
    %% compute the max norm error
    u_ex = u_exact(node);
    du = u_ex-u; ddu = 0*du;
    inpml = (x>-a + wpml + 10*eps).*(x<a -wpml - 10*eps)...
        .*(y>-a + wpml + 10*eps).*(y<a-wpml - 10*eps);
    idx0 = inpml>0;
    ddu(idx0) = du(idx0);
    du = du(idx0);
    rel_max_err(i) = norm(du,inf)/norm(u_ex,inf);
    max_err(i) = norm(du,inf);
    rel_l2_err(i) = norm(du,2)/norm(u_ex(idx0),2);
    l2_err(i) = norm(du,2)*h;
    toc;
end

figure(2);
subplot(2,2,1);
FJ_showrate(high_omega(1:n),max_err);
subplot(2,2,2);
FJ_showrate(high_omega(1:n),rel_max_err);
subplot(2,2,3);
FJ_showrate(high_omega(1:n),l2_err);
subplot(2,2,4);
FJ_showrate(high_omega(1:n),rel_l2_err);



       
%% Show result
% figure(3);
% FJ_showresult(node,elem,real(ddu))

fprintf('Max Error     ');  fprintf('  &  %.2e',rec_maxerr');  fprintf('\n\n');
fprintf('L2 Error      ');  fprintf('  &  %.2e',rec_l2err');  fprintf('\n\n');

    