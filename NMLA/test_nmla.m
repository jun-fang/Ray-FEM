% test the accuray order for one point source problem with NMLA corrected
% version

global omega xs ys;
xs = 2; ys = 2;

pde = Helmholtz_data_point_source;
Rest = 1;
Nray = 1;                  % one ray direction
sec_opt = 1;
plt = 0;



rec_omega = [40 80 160 320]*pi;
rec_error = rec_omega;
rec_r = rec_omega;
rec_cr = rec_r;

for rec_i = 1: length(rec_omega)
    
    omega = rec_omega(rec_i);
    NPW = 6;
    high_omega = omega;
    low_omega = sqrt(omega);
    h = 1/(NPW*round(high_omega/(2*pi)));
    ch = 1/(2*NPW*round(low_omega/(2*pi)));
    ds = 1/2;
    [node,elem] = squaremesh([-ds,ds,-ds,ds],h);
    [cnode,celem] = squaremesh([-ds,ds,-ds,ds],ch);
    ld = 2*ds;
    [lnode,lelem] = squaremesh([-ld,ld,-ld,ld],h);
    
    
    
    %% NMLA
    
    cN = size(cnode,1);
    cnumray_angle = zeros(cN,Nray);
    cnumray = zeros(cN,Nray);
    cr = zeros(cN,Nray);
    
    
    u_std = pde.ex_u(lnode);
    Du = pde.Du(lnode);
    ux = Du(:,1); uy = Du(:,2);
    
    fprintf('NMLA time: \n');
    tic;
    for i = 1:1 %cN
        x0 = cnode(i,1);  y0 = cnode(i,2);
        c0 = pde.speed(cnode(i,:));
        [cnumray_angle(i),r] = NMLA_2D_2nd(x0,y0,c0,omega,Rest,lnode,lelem,u_std,ux,uy,pde,1/5,Nray,'num',sec_opt,plt);
    end
    
    numray_angle = interpolation(cnode,celem,node,cnumray_angle);
    toc;
    
    exray_angle = ex_ray_angle(node,xs,ys);
    dang = numray_angle - exray_angle;
    error = norm(dang)*h;
    
    rec_error(rec_i) = error;
    rec_r(rec_i) = r;
    
    
    omega = low_omega;
    [cnumray_angle(i),r] = NMLA_2D_2nd(x0,y0,c0,omega,Rest,lnode,lelem,u_std,ux,uy,pde,1/5,Nray,'num',sec_opt,plt);
    rec_cr(rec_i) = r;
end

% showrate(rec_omega, rec_error)

