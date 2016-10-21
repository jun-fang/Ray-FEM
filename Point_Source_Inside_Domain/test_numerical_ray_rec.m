%% A test example for Ray-FEM: assembling with ray information

rec_omega = [ 50 100 200 400]*pi;
n = 4;
ray_maxerr = zeros(1,n);
ray_l2err = zeros(1,n);
nr_maxerr = zeros(1,n);
nr_l2err = zeros(1,n);
er_maxerr = zeros(1,n);
er_l2err = zeros(1,n);

for j = 1:n
    j
    tic;
    
    omega = rec_omega(j);
    xc = 0; yc = 0;
    aa = 1/2;
    sigma = 1/100;
    
    %% Load source data
    source = @(x) -4*sqrt(omega)*1i*1/(2*pi*sigma^2)*exp( -( (x(:,1)-xc).^2 + (x(:,2)-yc).^2  )/(2*sigma^2) );
    u_exact = @(p) sqrt(omega)*besselh(0,1,omega*sqrt((p(:,1)-xc).^2 + (p(:,2)-yc).^2));
    u_pl = @(x) (2/pi./(sqrt((x(:,1)-xc).^2 + (x(:,2)-yc).^2))).^0.5.*exp(1i*(omega*sqrt((x(:,1)-xc).^2 + (x(:,2)-yc).^2) - pi/4));
    speed = @(x) ones(size(x(:,1)));
    Nray = 1;
    
    %% Set up
    plt = 0;                   % show solution or not
    fquadorder = 3;            % numerical quadrature order
    NPW = 8;
    h = 1/round((NPW*omega)/(2*pi));
    
    wavelength = 2*pi/omega;
    wpml = 4*round(wavelength/h)*h;              % width of PML
    sigmaMax = 50/wpml;                % Maximun absorbtion
    
    a = aa + wpml;
    
    low_omega = 2*sqrt(omega);
    low_wl = 2*pi/low_omega;
    low_wpml = 4*ceil(low_wl/h)*h;
    low_sigmaMax = 50/low_wpml;
    
    %% Solving the low frequency equation
    low_a = a + ceil(low_wl/h)*h;
    [low_node,low_elem] = squaremesh([-low_a,low_a,-low_a,low_a],h);
    % [low_u] = Standard_FEM_PML_PointSource(low_node,low_elem,low_omega,low_wpml,low_sigmaMax,xc,yc,speed,fquadorder,plt);
    low_u = sqrt(low_omega)*besselh(0,1,low_omega*sqrt((low_node(:,1)-xc).^2 + (low_node(:,2)-yc).^2));
    % figure(1);
    % FJ_showresult(low_node,low_elem,real(low_u));
    
    %% Set up the high frequency computational domain
    
    bd = 1/omega^(1/4);
    bd = round(bd/(4*h))*4*h;
%     bd = 0.3;
    
    % a = 1/2+wpml;
    % [node,elem] = squaremesh([-a,a,-a,a],h);
    [node,elem] = squaremesh_annulus([-a,a,-a,a],[-bd,bd,-bd,bd],h);
    [cnode,celem] = squaremesh([-a,a,-a,a],4*h);
    % figure(2);
    % showmesh(cnode,celem);
    
    %% MPM to capture ray information
    cN = size(cnode,1);
    cnumang = zeros(cN,Nray);
    low_h = low_wl/NPW;
    tic;
    for i = 1:cN %round(N/100)
        x0 = cnode(i,1);  y0 = cnode(i,2);
        
        if (abs(x0)>=bd - 10*eps) || (abs(y0)>=bd - 10*eps)
            r = -low_wl:low_h:low_wl;  r = r';
            
            %% MPM sampling on one direction
            est_ang = pi/2;
            x = r*cos(est_ang) + x0;  y = r*sin(est_ang) + y0;
            xy = [x, y];
            u = interpolation(low_node,low_elem,xy,low_u);
            u = u.*sqrt((x-xc).^2 + (y-yc).^2).^(1/2);
            [z] = Matrix_Pencil(u,1);
            
            dang1 = real(acos(1i*imag(log(z))/(1i*low_omega*low_h)));
            ang1 = [est_ang + dang1, est_ang - dang1];
            ang1 = principal_angle(ang1);
            
            
            %% MPM sampling on one direction
            est_ang = pi/5;
            x = r*cos(est_ang) + x0;  y = r*sin(est_ang) + y0;
            xy = [x, y];
            u = interpolation(low_node,low_elem,xy,low_u);
            u = u.*sqrt((x-xc).^2 + (y-yc).^2).^(1/2);
            [z] = Matrix_Pencil(u,1);
            
            dang2 = real(acos(1i*imag(log(z))/(1i*low_omega*low_h)));
            ang2 = [est_ang + dang2, est_ang - dang2];
            ang2 = principal_angle(ang2);
            
            %% MPM sampling on another orthogonal direction
            est_ang =  3*pi/5;
            x = r*cos(est_ang) + x0;  y = r*sin(est_ang) + y0;
            xy = [x, y];
            u = interpolation(low_node,low_elem,xy,low_u);
            u = u.*sqrt((x-xc).^2 + (y-yc).^2).^(1/2);
            [z] = Matrix_Pencil(u,1);
            
            dang3 = real(acos(1i*imag(log(z))/(1i*low_omega*low_h)));
            ang3 = [est_ang + dang3, est_ang - dang3];
            ang3 = principal_angle(ang3);
            
            %% correction
            temp = find_angle_3(ang1, ang2, ang3, eps);
            if (length(temp) == 1)
                cnumang(i) = temp;
            else
                tol = pi/360;
                temp = find_angle_3(ang1, ang2, ang3, tol);
                while (length(temp) ~= 1)  && (tol < 5*pi/180)
                    tol = tol + pi/360;
                    temp = find_angle_3(ang1, ang2, ang3, tol);
                end
                cnumang(i) = mean(temp);
            end
        end
    end
    toc;
    
    cnumray = exp(1i*cnumang).*(1 - (abs(cnode(:,1))<bd-10*eps).*(abs(cnode(:,2))<bd-10*eps));
    cexray = ex_ray(cnode,xc,yc,1).*(1 - (abs(cnode(:,1))<bd-10*eps).*(abs(cnode(:,2))<bd-10*eps));
    cexang = ex_ray(cnode,xc,yc,0).*(1 - (abs(cnode(:,1))<bd-10*eps).*(abs(cnode(:,2))<bd-10*eps));
    
    numray = interpolation(cnode,celem,node,cnumray);
    numang = atan2(imag(numray), real(numray));
    numray = exp(1i*numang);
    
    %% Exact ray information
    xx = node(:,1)-xc; yy = node(:,2)-yc;
    ray = atan2(yy,xx);
    ray = exp(1i*ray);
    exray = ray.*(1 - (abs(xx)<eps).*(abs(yy)<eps));
    exrayang = ex_ray(node,xc,yc,0);
    
    dray = exray - numray;
    dang = exrayang - numang;
    ray_maxerr(j) = norm(dray,inf);
    ray_l2err(j) = norm(dray,2)/norm(exray,2);
    
    N = size(node,1);       % number of grid points
    Nray = size(ray,2);     % number of rays crossing at each grid node
    Ndof = N*Nray;
    
    
    %% Assembling the linear system
    ray = exray;
    eA = assemble_Helmholtz_matrix_with_ray(node,elem,omega,wpml,sigmaMax,speed,ray,fquadorder);
    eb = assemble_RHS_with_ray(node,elem,omega,source,speed,ray,fquadorder);
    
    ray = numray;
    nA = assemble_Helmholtz_matrix_with_ray(node,elem,omega,wpml,sigmaMax,speed,ray,fquadorder);
    nb = assemble_RHS_with_ray(node,elem,omega,source,speed,ray,fquadorder);
    
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
    
    
    %% Numerical Ray-FEM
    ray = numray;    
    dx = real(ray(idx2)).*BdNode2(:,1) + imag(ray(idx2)).*BdNode2(:,2);
    v = zeros(Ndof,1);
    v(idx2) = u_exact(BdNode2).*exp(-1i*omega*dx);
    nb = nb - nA*v;
    v(freeNode) = nA(freeNode,freeNode)\nb(freeNode);
    
    grad = ray(:);
    grad = [real(grad),imag(grad)];
    repnode = repmat(node,Nray,1);
    temp = grad(:,1).*repnode(:,1) + grad(:,2).*repnode(:,2);
    
    c = speed(node);    % medium speed
    k = omega./c;           % wavenumber
    kk = repmat(k,1,Nray);
    
    nu = v.*exp(1i*kk(:).*temp);
    nu = reshape(nu,N,Nray);
    nu = sum(nu,2);
    
    
    %% Exact Ray-FEM
    ray = exray;    
    dx = real(ray(idx2)).*BdNode2(:,1) + imag(ray(idx2)).*BdNode2(:,2);
    v = zeros(Ndof,1);
    v(idx2) = u_exact(BdNode2).*exp(-1i*omega*dx);
    eb = eb - eA*v;
    v(freeNode) = eA(freeNode,freeNode)\eb(freeNode);
    
    grad = ray(:);
    grad = [real(grad),imag(grad)];
    repnode = repmat(node,Nray,1);
    temp = grad(:,1).*repnode(:,1) + grad(:,2).*repnode(:,2);
    
    c = speed(node);    % medium speed
    k = omega./c;           % wavenumber
    kk = repmat(k,1,Nray);
    
    eu = v.*exp(1i*kk(:).*temp);
    eu = reshape(eu,N,Nray);
    eu = sum(eu,2);
    
    
    %% compute the max norm error
    u_ex = u_exact(node);
    inpml = (x>-a + wpml + 10*eps).*(x<a -wpml - 10*eps)...
        .*(y>-a + wpml + 10*eps).*(y<a-wpml - 10*eps);
    idx0 = inpml>0;
    ndu = u_ex-nu; edu = u_ex-eu;
    ndu = ndu(idx0); edu = edu(idx0);
    nr_maxerr(j) = norm(ndu,inf);
    nr_l2err(j) = norm(ndu,2)/norm(u_ex(idx0),2);
    er_maxerr(j) = norm(edu,inf);
    er_l2err(j) = norm(edu,2)/norm(u_ex(idx0),2);
    toc;
end

fprintf('Ray    Max Error     ');  fprintf('  &  %.2e',ray_maxerr');  fprintf('\n\n');
fprintf('Ray    L2 Error      ');  fprintf('  &  %.2e',ray_l2err');  fprintf('\n\n');

fprintf('NR-FEM Max Error     ');  fprintf('  &  %.2e',nr_maxerr');  fprintf('\n\n');
fprintf('NR-FEM L2 Error      ');  fprintf('  &  %.2e',nr_l2err');  fprintf('\n\n');

fprintf('ER-FEM Max Error     ');  fprintf('  &  %.2e',er_maxerr');  fprintf('\n\n');
fprintf('ER-FEM L2 Error      ');  fprintf('  &  %.2e',er_l2err');  fprintf('\n\n');






%% Show result
% figure(2);
% FJ_showresult(node,elem,real(u))
% 
% 
% figure(3);
% FJ_showresult(node,elem,real(ddu))
