% test the accuray order for one point source problem in homogeneous medium

xs = 2; ys = 2;  
speed = @(p) ones(size(p(:,1)));

Rest = 1;
Nray = 1;                  % one ray direction
sec_opt = 0;
plt = 0;
NPW = 6;
ds = 1/2;
ld = 2*ds;

omegas = [ 80 160 320 640]*pi;
% omega1 = omegas + 50;
% omega2 = omegas - 50;
omega1 = 1.5*omegas;
omega2 = 0.75*omegas;
errs = omegas;
cerrs = omegas;
cerrs1 = omegas;
cerrs2 = omegas;



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
    
    %% case 1
    omega = omega1(oi);
    ul = besselh(0,1,omega*rr);
    [ux,uy] = num_derivative(ul,h,2);
    
    fprintf('NMLA time: \n');
    tic;
    for i = 1:cN
        x0 = cnode(i,1);  y0 = cnode(i,2);
        c0 = speed(cnode(i,:));
        [cnumray_angle(i),r] = NMLA(x0,y0,c0,omega,Rest,lnode,lelem,ul,ux,uy,[],1/5,Nray,'num',sec_opt,plt);
    end
    cexray_angle = ex_ray(cnode,xs,ys,0);
    cerrs1(oi) = norm(cnumray_angle - cexray_angle,inf);
    ang1 = cnumray_angle;
    
    %% case 2
    omega = omega2(oi);
    ul = besselh(0,1,omega*rr);
    [ux,uy] = num_derivative(ul,h,2);
    
    fprintf('NMLA time: \n');
    tic;
    for i = 1:cN
        x0 = cnode(i,1);  y0 = cnode(i,2);
        c0 = speed(cnode(i,:));
        [cnumray_angle(i),r] = NMLA(x0,y0,c0,omega,Rest,lnode,lelem,ul,ux,uy,[],1/5,Nray,'num',sec_opt,plt);
    end
    cexray_angle = ex_ray(cnode,xs,ys,0);
    cerrs2(oi) = norm(cnumray_angle - cexray_angle,inf);
    ang2 = cnumray_angle;
    
    ang = (sqrt(omega1(oi)).*ang1 - sqrt(omega2(oi)).*ang2)...
        ./(sqrt(omega1(oi)) - sqrt(omega2(oi)));
    
%     ang = (omega1(oi).*ang1 - omega2(oi).*ang2)...
%         ./(omega1(oi) - omega2(oi));
    
    cnumray_angle = ang;
    
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


figure(3);
subplot(1,2,1);
show_convergence_rate(omegas,cerrs1,'omega');
subplot(1,2,2);
show_convergence_rate(omegas,cerrs2,'omega');

