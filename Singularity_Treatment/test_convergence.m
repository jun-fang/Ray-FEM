%% Test the convergence order of SFEM and RayFEM for singularity removel problem


%% S-FEM h convergence:   error \sim O(h^2)
if (0)
    xs = 0; ys = 0;
    omega = 40*pi;
    wl = 2*pi/omega;
    epsilon = 0.143;
    speed = @(p) ones(size(p(:,1)));
    
    wpml = 0.2;
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



%% SFEM omega convergence: error \sim O(\omega^2.5)
if (0)
    xs = 0; ys = 0;
    omega = 40*pi;
    wl = 2*pi/omega;
    epsilon = 0.143;
    speed = @(p) ones(size(p(:,1)));
    
    wpml = 0.1;
    sigmaMax = 25/wpml;
    fquadorder = 3;
    a = 1/2;
    
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
    
    
    figure(12);
    subplot(1,2,1);
    show_convergence_rate(omegas,rhss,'omega',[],'||f||_{L^2(\Omega)}');
    subplot(1,2,2);
    show_convergence_rate(omegas,errors,'omega',[],'||u - u_h||_{L^2(\Omega)}');
    
end






%% Ray-FEM h convergence:  error \sim O(h^2)
if(0)
    xs = 0; ys = 0;
    omega = 40*pi;
    wl = 2*pi/omega;
    epsilon = 0.143;
    speed = @(p) ones(size(p(:,1)));
    
    wpml = 0.2;
    sigmaMax = 25/wpml;
    fquadorder = 3;
    a = 1/2;
    
    nt = 4;
    hs = zeros(nt,1);
    RFEM_errs = zeros(nt,1);  % error on physical domain
    h = 1/100;
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
        
        [u,~,~,v] = RayFEM_singularity(node,elem,xs,ys,omega,epsilon,wpml,sigmaMax,ray,speed,@sing_rhs_homo,fquadorder);
        
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
        
        [err, rel_L2_err] = RayFEM_smooth_solution_error(node,elem,xs,ys,omega,epsilon,wpml,ray,speed,v,9);
        RFEM_errs(ii) = err;%norm(du_phy)*h;
        toc;
    end
    
    figure(21);
    show_convergence_rate(hs,RFEM_errs,'h',[],'||u - u_h||_{L^2(\Omega)}');
end



%% Ray FEM omega convergence: error \sim O(\omega^2.5)???
if (1)
    xs = 0; ys = 0;
%     omega = 40*pi;
    
    epsilon = 0.143;
    speed = @(p) ones(size(p(:,1)));
    
    wpml = 0.1;
    sigmaMax = 25/wpml;
    fquadorder = 3;
    a = 1/2;
    
    
    nt = 5;
    errors = zeros(1,nt);
    rhss = zeros(1,nt);
    omegas = pi*[160,240,320,480,640];
    NPW = 4;
    
    
    
    for ii = 1:nt
        ii
        omega = omegas(ii);
        wl = 2*pi/omega;
        h = 1/round(1/(wl/NPW));
        1/h
        epsilon = sqrt(18/omega);
        epsilon
    
%     h = 1/240;
    
    [node,elem] = squaremesh([-a,a,-a,a],h);
    
    %% Exact ray information
    xx = node(:,1)-xs;  yy = node(:,2)-ys;
    rr = sqrt(xx.^2 + yy.^2);
    ray = atan2(yy,xx);
    ray = exp(1i*ray).*(rr>10*eps);
    
    
    
%         omega = 40*pi;
        rhs = sing_rhs_homo(epsilon,omega,node,xs,ys);
        rhss(ii) = norm(rhs)*h;
        
     
        
        
        tic;
        [u,~,~,v] = RayFEM_singularity(node,elem,xs,ys,omega,epsilon,wpml,sigmaMax,ray,speed,@sing_rhs_homo,fquadorder);
        
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
        
        dd = 0*du; dd(idx) = du(idx);
%         figure(2);showsolution(node,elem,real(dd));
        
        
        [err, rel_L2_err] = RayFEM_smooth_solution_error(node,elem,xs,ys,omega,epsilon,wpml,ray,speed,v,9);
        
        errors(ii) = err;%norm(du_phy)*h;%/norm(uex(idx));
        toc;
    end
    
    
    figure(22);
    subplot(1,2,1);
    show_convergence_rate(omegas(1:nt),rhss,'omega',[],'||f||_{L^2(\Omega)}');
    subplot(1,2,2);
    show_convergence_rate(omegas(1:nt),errors,'omega',[],'||u - u_h||_{L^2(\Omega)}');
    
end



