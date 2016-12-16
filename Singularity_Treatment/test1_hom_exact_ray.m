xs = 0; ys = 0;

epsilon = 50/(80*pi);
speed = @(p) ones(size(p(:,1)));

wpml = 0.1;
sigmaMax = 25/wpml;
fquadorder = 3;
a = 1/2;


nt = 3;
errors = zeros(1,nt);
rhss = zeros(1,nt);
omegas = pi*[120,160,240,320,480,640];
NPW = 4;



for ii = 1:nt
    ii
    omega = omegas(ii);
    wl = 2*pi/omega;
    h = 1/round(1/(wl/NPW));
    1/h
%     epsilon = sqrt(18/omega);
    
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
