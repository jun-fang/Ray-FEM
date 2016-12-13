% test the accuray order for one point source problem in homogeneous medium

xs = 2; ys = 2;  
speed = @(p) ones(size(p(:,1)));

pde = Helmholtz_data1_0;
Rest = 1;
Nray = 1;                  % one ray direction
sec_opt = 0;
plt = 0;
NPW = 6;
ds = 1/2;
ld = 2*ds;

omegas = [40 80 160 320 640]*pi;
errs = omegas;
cerrs = omegas;

for oi = 1: length(omegas)
    
    omega = omegas(oi);
    h = 1/(NPW*round(omega/(2*pi)));
    ch = 1/(10*round(NPW*sqrt(omega)/(2*pi)/4));
    
    [node,elem] = squaremesh([-ds,ds,-ds,ds],h);
    [cnode,celem] = squaremesh([-ds,ds,-ds,ds],ch);
    [lnode,lelem] = squaremesh([-ld,ld,-ld,ld],h);
  
    
    %% NMLA
    cN = size(cnode,1);
    cnumray_angle = zeros(cN,Nray);
    
    x = lnode(:,1); y = lnode(:,2);
    rr = sqrt((x-xs).^2 + (y-ys).^2);
    ul = besselh(0,1,omega*rr);
    [ux,uy] = num_derivative(ul,h,2);
    
    fprintf('NMLA time: \n');
    tic;
    for i = 1:cN
        x0 = cnode(i,1);  y0 = cnode(i,2);
        c0 = speed(cnode(i,:));
        [cnumray_angle(i),r] = NMLA(x0,y0,c0,omega,Rest,lnode,lelem,ul,ux,uy,pde,1/5,Nray,'num',sec_opt,plt);
    end
    cexray_angle = ex_ray(cnode,xs,ys,0);
    cerrs(oi) = norm(cnumray_angle - cexray_angle,inf);

    numray_angle = interpolation(cnode,celem,node,cnumray_angle);
    toc;
    
    exray_angle = ex_ray(node,xs,ys,0);
    errs(oi) = norm(numray_angle - exray_angle,inf);
    
end

figure(2);
subplot(1,2,1);
show_convergence_rate(omegas,errs,'omega');
subplot(1,2,2);
show_convergence_rate(omegas,cerrs,'omega');

