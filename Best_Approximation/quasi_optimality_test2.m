global omega a xs ys;
pde = Helmholtz_data2;
fquadorder = 9;
sec_opt = 0;
plt = 0;
Nray = 4;
NPW = 6;
a = 1/2; xs = 2; ys = 2;
Rest = 2.5*sqrt(2);

ome = 0:2;
omegas = pi*(20+20*ome); 
ray_err = 0*omegas;
solu_err = 0*omegas;
best_err = 0*omegas;

for ii = 1:length(omegas)
    
    omega = omegas(ii);
    fprintf(['\n' '-'*ones(1,80) '\n']);
    fprintf('Case: %d, omega/2pi = %d \n',ii,round(omega/(2*pi)));
    
    h = 1/(NPW*round(omega/(2*pi)));
    ch = 1/(10*round(1.5*sqrt(omega)/(2*pi)));
    
    [node,elem] = squaremesh([-a,a,-a,a],h);
    [cnode,celem] = squaremesh([-a,a,-a,a],ch);
    
    cN = size(cnode,1);
    cnumray = zeros(cN,Nray);
    u =[];ux=[];uy=[];
    
    fprintf('NMLA time: \n');
    tic;
    for i = 1:cN
        x0 = cnode(i,1);  y0 = cnode(i,2);
        c0 = pde.speed(cnode(i,:));
        Rest = 100 ;
        [cnumray(i,:)] = NMLA_2D_correction(x0,y0,c0,omega,Rest,node,elem,u,ux,uy,pde,1/8,Nray,'ex',sec_opt,plt);
    end
%     cnumray = cnumray + 0.02*(0.5-rand(size(cnumray)));
    cnumray = exp(1i*cnumray);
    numray = interpolation(cnode,celem,node,cnumray);
    toc;
    
    exray = pde.ray(node);
    dray = exray-numray;
    ray_err(ii) = norm(dray(:))*h;
    ray = numray;
    
%     fprintf('Ray-FEM time: \n');
%     tic;
%     [uh] = Ray_FEM_IBC_1(node,elem,omega,pde,ray,fquadorder,plt);
%     ue = pde.ex_u(node);
%     solu_err(ii) = norm(uh-ue)*h;
%     toc;

    fprintf('Best approximation time: \n');
    tic;
    [error] = best_approximation(node,elem,ray,omega,pde);
    best_err(ii) = error;
    toc;
end

% figure(1);
% plot(omegas,solu_err./best_err,'*-');
