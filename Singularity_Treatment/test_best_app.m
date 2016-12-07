xs = 0; ys = 0;

epsilon = sqrt(18/40/pi);
speed = @(p) ones(size(p(:,1)));


a = 1;


nt = 5;
errors1 = zeros(1,nt);
errors2 = zeros(1,nt);
omegas = pi*[40,60,80,120,160];
NPW = 6;


for ii = 1:nt
    ii
    tic;
    omega = omegas(ii);
    wl = 2*pi/omega;
    h = 1/round(1/(wl/NPW));
    1/h
        
    [node,elem] = squaremesh([-a,a,-a,a],h);
    
    %% Exact ray information
    xx = node(:,1)-xs;  yy = node(:,2)-ys;
    rr = sqrt(xx.^2 + yy.^2);
    ray = atan2(yy,xx);
    ray = exp(1i*ray).*(rr>10*eps);
    
    grad = ray(:);
    grad = [real(grad),imag(grad)];
    temp = grad(:,1).*node(:,1) + grad(:,2).*node(:,2);
    
    k = omega./speed(node);           % wavenumber
    
    ub = 1i/4*besselh(0,1,omega*rr);
    cf = cutoff(epsilon,2*epsilon,node,xs,ys);
    uex = (1-cf).*ub;
    uex(rr<epsilon) = 0;
    
    v = uex.*exp(-1i*k.*temp);
    
    
    
    %% high order quadrature
    fquadorder = 9;
    %% FEM set up
    N = size(node,1);       % number of grid points
    NT = size(elem,1);      % number of triangle elements
    Nray = size(ray,2);     % number of rays crossing at each grid node
    
    c = speed(node);    % medium speed
    k = omega./c;           % wavenumber
    kk = repmat(k,1,Nray);
    
    xmax = max(node(:,1)); xmin = min(node(:,1));
    ymax = max(node(:,2)); ymin = min(node(:,2));
    
    
    %% Numerical Quadrature
    [lambda,weight] = quadpts(fquadorder);
    phi = lambda;           % the value of piesewise linear basis functions at numerical quadrature nodes
    nQuad = size(lambda,1);
    
    
    %% Compute geometric quantities and gradient of local basis
    [~,area] = gradbasis(node,elem);
    
    
    %% Get relative L2 error
    err = zeros(NT,1);
    rel_err = err;
    exp_phaseii = zeros(NT,3,Nray);
    for p = 1:nQuad
        pxy = lambda(p,1)*node(elem(:,1),:) ...
            + lambda(p,2)*node(elem(:,2),:) ...
            + lambda(p,3)*node(elem(:,3),:);
        uhp = zeros(NT,1);
        for ni = 1: Nray
            for i = 1:3
                ray_ni = ray(elem(:,i),ni);
                gradtempi = [real(ray_ni), imag(ray_ni)];
                fphasei = gradtempi(:,1).*pxy(:,1) + gradtempi(:,2).*pxy(:,2);
                ki = kk(elem(:,i),ni);
                exp_phaseii(:,i,ni) = exp(1i*ki.*fphasei);
            end
            nii = (ni - 1)*N;
            uhp = uhp + v(elem(:,1) + nii)*phi(p,1).*exp_phaseii(:,1,ni) + ...
                v(elem(:,2) + nii)*phi(p,2).*exp_phaseii(:,2,ni) + ...
                v(elem(:,3) + nii)*phi(p,3).*exp_phaseii(:,3,ni);
        end
        
        x = pxy(:,1); y = pxy(:,2);
        rr = sqrt((x-xs).^2 + (y-ys).^2);
        
        ub = 1i/4*besselh(0,1,omega*rr);
        cf = cutoff(epsilon,2*epsilon,pxy,xs,ys);
        uex = (1-cf).*ub;
        uex(rr<epsilon) = 0;
        
        
        abs_err = abs(uex - uhp);
        err = err + weight(p)*(abs_err).^2;
        rel_err = rel_err + weight(p)*(abs(uex)).^2;
    end
    
    err = area.*err;
    rel_err = area.*rel_err;
    err(isnan(err)) = 0;
    err = sqrt(sum(err));
    rel_err = sqrt(sum(rel_err));
    rel_L2_err = err./rel_err;
    errors1(ii) = err;
    errors2(ii) = rel_L2_err;
    toc;
    
end


figure(10);
subplot(1,2,1);
show_convergence_rate(omegas(1:nt),errors1);
subplot(1,2,2);
show_convergence_rate(omegas(1:nt),errors2);