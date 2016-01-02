global omega;
pde = Helmholtz_data8;
a = 1/2;
plt = 0;
Nray = 1;
Rest = 1;
omega = 60*pi;
% omega = sqrt(omega)

xs = 1/8;
ys = 1/10;

NPW = 20;
h = 1/(NPW*round(omega/(2*pi)));
a = 1;
[node,elem] = squaremesh([-a,a,-a,a],h);
ch = 1/60;
a = 1/2;
[cnode,celem] = squaremesh([-a,a,-a,a],ch);
    cN = size(cnode,1);
    cnumray = zeros(cN,Nray);

    cray = ex_ray(cnode,xs,ys);

    u1 = pde.ex_u(node);
    Du = pde.Du(node);
    ux = Du(:,1);
    uy = Du(:,2);
    
    for i = 1:cN
    x0 = cnode(i,1);  y0 = cnode(i,2);
    d0 = sqrt((x0-xs)^2 + (y0-ys)^2);
    if d0 <= r
        cnumray(i,:) = cray(i,:);
    else
        Rest = d0;
        c0 = pde.speed(cnode(i,:));
        [cnumray(i,:)] = NMLA_2D_2nd(x0,y0,c0,omega,Rest,node,elem,u1,ux,uy,pde,pct,Nray,'num',opt,plt);
    end
    end
    


    
    
    
    cdiffang = angle_error(cnumray,cray);
norm(cdiffang,2)/norm(cray)
norm(cdiffang,inf)






if (0)
clear;
pde = Helmholtz_data7;
fquadorder = 6;
plt = 0;
Nray = 1;
global omega;
global a;
lg_a = 0.6;
sm_a = 1/2;
n = 1;
dfang = zeros(1,n);
dcang = dfang;
ss = dfang;
data = 'num';
Rest = 2.5;
omega = 40*pi;

NPW = 10;

 cp_omega = [10 20 30 40]*2*pi;

for i = 1:n
    i
%     omega = cp_omega(i);
%     omega = 2*omega;
    ss(i) = omega;
    
    high_omega = omega;
    low_omega = sqrt(high_omega);
    
    %     omega = 2*low_omega;
    omega = low_omega;
    
    h = 1/(NPW*round(high_omega/(2*pi)));
    ch = 1/(NPW*round((low_omega)/(2*pi)));
    
    a = lg_a;
    [lnode,lelem] = squaremesh([-a,a,-a,a],h);
    lN = size(lnode,1);
    
    u = pde.ex_u(lnode);
    [ux,uy] = num_derivative(u,h,4);
%     du = pde.Du(lnode);
%     rand1 = rand(size(u));
%     rand2 = rand(size(du));
%     u = u + rand1*norm(u,inf)*10^(-2);
%     du = du + rand2*norm(du,inf)*10^(-2);
%     ux = du(:,1);
%     uy = du(:,2);
    
    
    a = sm_a;
    [node,elem] = squaremesh([-a,a,-a,a],h);
    ray_dir = pde.ray_ang(node);
    
    [cnode,celem] = squaremesh([-a,a,-a,a],ch);
    cN = size(cnode,1);
    
    cnumray = zeros(cN,Nray);
    
    cray_dir = pde.ray_ang(cnode);
    
    
    for j = 1: 1 %cN
        x0 = cnode(j,1);  y0 = cnode(j,2);
        c0 = pde.speed(cnode(j,:));
%         x0 = 0; y0 = 0; c0 = 1;
        [cnumray(j,:)] = NMLA_2D_2nd(x0,y0,c0,omega,Rest,lnode,lelem,u,ux,uy,pde,1/5,Nray,data,plt);
    end
    
    cnumray(1,:)/pi - [0.25 0.75 1.25 1.75]
    
    diffcang = cnumray -cray_dir;
    dcang(i) = ch*norm(diffcang(:),2)/(ch*norm(cray_dir(:),2));
    numray = interpolation(cnode,celem,node,cnumray);
    diffang = numray - ray_dir;
    dfang(i) = h*norm(diffang(:),2)/(h*norm(ray_dir(:),2));
    norm(diffcang(:),inf)
    norm(diffang(:),inf)
end

% numray = exp(1i*numray);
% [~,~,~,~,rel_L2_err1] = Ray_FEM_IBC_1(node,elem,omega,pde,numray,fquadorder,plt);
% % 
ray = pde.ray(node);
[~,A,v,b,rel_L2_err2] = Ray_FEM_IBC_1(node,elem,omega,pde,ray,fquadorder,plt);

% numray = ray_dir + rand(size(ray_dir))*0.0055;
% numray = exp(1i*numray);
% [~,~,~,~,rel_L2_err3] = Ray_FEM_IBC_1(node,elem,omega,pde,numray,fquadorder,plt)


% ss
% dcang
% dfang
% showrate(ss,dfang)
end