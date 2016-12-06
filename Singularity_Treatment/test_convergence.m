%% Test the convergence order of SFEM and RayFEM for singularity removel problem




%% SFEM: omega convergence
if (0)
xs = 0; ys = 0;
omega = 40*pi;
wl = 2*pi/omega;
epsilon = sqrt(18/40/pi);
speed = @(p) ones(size(p(:,1)));

wpml = 0.1;
sigmaMax = 25/wpml;
fquadorder = 3;
a = 1;

h = 1/1000;

[node,elem] = squaremesh([-a,a,-a,a],h);

%% Exact ray information
xx = node(:,1)-xs;  yy = node(:,2)-ys;
rr = sqrt(xx.^2 + yy.^2);
ray = atan2(yy,xx);
ray = exp(1i*ray).*(rr>10*eps);

nt = 4;
errors = zeros(1,nt);
rhss = zeros(1,nt);
omegas = pi*[40,60,80,100];
for ii = 1:4
    ii
    omega = omegas(ii);
    
    rhs = sing_rhs_homo(epsilon,omega,node,xs,ys);
    rhss(ii) = norm(rhs)*h;

    tic;
    % [~,A,b,~] = RayFEM_singularity(node,elem,xs,ys,omega,epsilon,wpml,sigmaMax,ray,speed,@sing_rhs_homo,fquadorder);
    [u0,A0,b0,~] = RayFEM_singularity(node,elem,xs,ys,omega,epsilon,wpml,sigmaMax,0*ray,speed,@sing_rhs_homo,fquadorder);
    
    x = node(:,1); y = node(:,2);
    rr = sqrt((x-xs).^2 + (y-ys).^2);
    
    ub = 1i/4*besselh(0,1,omega*rr);
    cf = cutoff(epsilon,2*epsilon,node,xs,ys);
    uex = (1-cf).*ub;
    uex(rr<epsilon) = 0;
    du = u0 - uex;
    
    idx = find( (x<=max(x)-wpml).*(x>= min(x)+wpml)...
        .*(y<= max(y)-wpml).*(y>= min(y)+wpml) );
    du_phy = du(idx);
    errors(ii) = norm(du_phy)*h;%/norm(uex(idx));
    toc;
end


figure(1);
subplot(1,2,1);
show_convergence_rate(omegas,rhss,'omega',[],'||f||_{L^2(\Omega)}');
subplot(1,2,2);
show_convergence_rate(omegas,erros,'omega',[],'||u - u_h||_{L^2(\Omega)}');

end






%% test h convergence: Ray-FEM error \sim O(h^2)
if(0)
    xs = 0; ys = 0;
    omega = 40*pi;
    wl = 2*pi/omega;
    epsilon = 0.143;
    speed = @(p) ones(size(p(:,1)));
    
    wpml = 25/128;
    sigmaMax = 25/wpml;
    fquadorder = 3;
    a = 1/2;
    
    nt = 3;
    hs = zeros(nt,1);
    RFEM_errs = zeros(nt,1);  % error on physical domain
    h = 1/200;
    for ii = 1:nt
        ii
        tic;
        h = h/2;
        hs(ii) = h;
        [node,elem] = squaremesh([-a,a,-a,a],h);
        
        %% Exact ray information
        xx = node(:,1)-xs;  yy = node(:,2)-ys;
        rr = sqrt(xx.^2 + yy.^2);
        ray = atan2(yy,xx);
        ray = exp(1i*ray).*(rr>10*eps);
        
        A = assemble_Helmholtz_matrix_with_ray_1(node,elem,omega,wpml,sigmaMax,speed,ray,fquadorder);
        b = assemble_RHS_with_sing_RayFEM(node,elem,xs,ys,omega,epsilon,wpml,sigmaMax,ray,speed,@sing_rhs_homo,fquadorder);
        
        %% Boundary conditions
        [bdNode,~,isBdNode] = findboundary(elem);
        freeNode = find(~isBdNode);
        N = size(node,1);  Nray = size(ray,2); n = round(sqrt(N));
        v = zeros(N,1);
        v(freeNode) = A(freeNode,freeNode)\b(freeNode);
        
        %% construct solution
        grad = ray(:);
        grad = [real(grad),imag(grad)];
        repnode = repmat(node,Nray,1);
        temp = grad(:,1).*repnode(:,1) + grad(:,2).*repnode(:,2);
        
        k = omega./speed(node);           % wavenumber
        kk = repmat(k,1,Nray);
        u = v.*exp(1i*kk(:).*temp);
        u = reshape(u,N,Nray);
        u = sum(u,2);
        
        x = node(:,1); y = node(:,2);
        rr = sqrt((x-xs).^2 + (y-ys).^2);
        
        ub = 1i/4*besselh(0,1,omega*rr);
        cf = cutoff(epsilon,2*epsilon,node,xs,ys);
        uex = (1-cf).*ub;
        uex(rr<epsilon) = 0;
        du = u - uex;
        
        idx = find( (x<=max(x)-wpml).*(x>= min(x)+wpml)...
            .*(y<= max(y)-wpml).*(y>= min(y)+wpml) );
        du_phy = du(idx);
        RFEM_errs(ii) = norm(du_phy)/norm(uex(idx));
        toc;
    end
    
    figure(10);
    show_convergence_rate(hs,RFEM_errs,'h',,[],'||u - u_h||_{L^2(\Omega)}');
end



%% test h convergence: S-FEM  error \sim O(h^2)
if (0)
    xs = 0; ys = 0;
    omega = 40*pi;
    wl = 2*pi/omega;
    epsilon = 0.143;
    speed = @(p) ones(size(p(:,1)));
    
    wpml = 25/128;
    sigmaMax = 25/wpml;
    fquadorder = 3;
    a = 1/2;
    
    nt = 4;
    hs = zeros(nt,1);
    SFEM_errs = zeros(nt,1);  % error on physical domain
    h = 1/100;
    for ii = 1:nt
        ii
        tic;
        h = h/2;
        hs(ii) = h;
        [node,elem] = squaremesh([-a,a,-a,a],h);
        A = assemble_Helmholtz_matrix_SFEM(node,elem,omega,wpml,sigmaMax,speed,fquadorder);
        b = assemble_RHS_with_sing_SFEM(node,elem,xs,ys,omega,wpml,sigmaMax,epsilon,fquadorder);
        
        %% Boundary conditions
        [bdNode,~,isBdNode] = findboundary(elem);
        freeNode = find(~isBdNode);
        N = size(node,1);  n = round(sqrt(N));
        u = zeros(N,1);
        u(freeNode) = A(freeNode,freeNode)\b(freeNode);
        
        x = node(:,1); y = node(:,2);
        rr = sqrt((x-xs).^2 + (y-ys).^2);
        
        ub = 1i/4*besselh(0,1,omega*rr);
        cf = cutoff(epsilon,2*epsilon,node,xs,ys);
        uex = (1-cf).*ub;
        uex(rr<epsilon) = 0;
        du = u - uex;
        
        idx = find( (x<=max(x)-wpml).*(x>= min(x)+wpml)...
            .*(y<= max(y)-wpml).*(y>= min(y)+wpml) );
        du_phy = du(idx);
        SFEM_errs(ii) = norm(du_phy)/norm(uex(idx));
        toc;
    end
    
    figure(11);
    show_convergence_rate(hs,SFEM_errs,'h',[],'||u - u_h||_{L^2(\Omega)}');
end




