global omega a;
pde = Helmholtz_data1;
fquadorder = 3;
sec_opt = 1;
plt = 0;
Nray = 1;
NPW = 6;
a = 1/2;
Rest = 2.5*sqrt(2);

ome = 0:30;
omegas = pi*(20+10*ome); 
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
        Rest = sqrt((x0-2)^2 + (y0-2)^2);
        [cnumray(i)] = NMLA_2D_2nd(x0,y0,c0,omega,Rest,node,elem,u,ux,uy,pde,1/5,Nray,'ex',sec_opt,plt);
    end
    cnumray = exp(1i*cnumray);
    numray = interpolation(cnode,celem,node,cnumray);
    toc;
    
    exray = ex_ray(node,2,2,1);
    ray_err(ii) = norm(exray-numray)*h;
    
    fprintf('Ray-FEM time: \n');
    tic;
    [uh] = Ray_FEM_IBC_1(node,elem,omega,pde,numray,fquadorder,plt);
    ue = pde.ex_u(node);
    solu_err(ii) = norm(uh-ue)*h;
    toc;

    fprintf('Best approximation time: \n');
    tic;
    [error] = best_approximation(node,elem,numray,omega,pde);
    best_err(ii) = error;
    toc;
end


save('qoc_1.mat','omegas', 'ray_err', 'solu_err', 'best_err');

figure(1);
plot(omegas,solu_err./best_err,'*-');

% qoc = [2.14383853061539,2.12906260837922,2.08161260494130,2.08084225786493,2.12558780001491,2.12255228939347,2.12230660788605,2.12429096852222,2.12426947685402,2.13323187551869,2.13077759960506,2.12930744773401,2.12636569770460,2.09979352655347,2.12092590632039,2.10316219597566]
