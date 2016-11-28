global omega a xs ys;
pde = Helmholtz_data2;
omega = 0*pi;
a = 1/2; xs = 2; ys = 2;
fquadorder= 9; plt = 0; NPW = 4; sec_opt = 0;
k = 2;
test_n = 8;
ray_err = zeros(1,test_n);
solu_err = ray_err;
best_err = ray_err;

omegas = pi*[10,20,30,40,50,60,70,80];

for oi = 1:test_n
    omega = omegas(oi);% + 10*pi;
    fprintf(['\n' '-'*ones(1,80) '\n']);
    fprintf('Case: %d, omega/2pi = %d \n',oi,round(omega/(2*pi)));
    
    h = 1/(NPW*round(omega/(2*pi)));
    [node,elem] = squaremesh([-a,a,-a,a],h);
    N = size(node,1); NT = size(elem,1); area = h*h/2;
       
    ch = 1/(10*round(1.5*sqrt(omega)/(2*pi)));
    [cnode,celem] = squaremesh([-a,a,-a,a],ch);
    
    cN = size(cnode,1);
    cnumray = zeros(cN,Nray);
    u =[];ux=[];uy=[];
    
    fprintf('NMLA time: \n');
    tic;
    for i = 1:cN
        x0 = cnode(i,1);  y0 = cnode(i,2);
        c0 = pde.speed(cnode(i,:));
        Rest = 100-sqrt(omega)/2;
        [cnumray(i,:)] = NMLA_2D_correction(x0,y0,c0,omega,Rest,node,elem,u,ux,uy,pde,1/8,Nray,'ex',sec_opt,plt);
    end
    cnumray = exp(1i*cnumray);
    numray = interpolation(cnode,celem,node,cnumray);
    toc;
    
    exray = pde.ray(node);
    dray = exray-numray;
    ray_err(oi) = norm(dray(:))*h;
    
    
    ray = numray;
%     ray = exray;
    
      
    %% Best approximation 
    if (1)
    fprintf('Best_approximation time: \n');
    tic;
    rnode = node; relem = elem; fh = h;
    for i = 1:k
        [rnode,relem] = uniformrefine(rnode,relem);
        fh = fh/2;
    end
    
    xnode = rnode(:,1); ynode = rnode(:,2);
    rN = size(rnode,1);
    Nray = size(ray,2);
    ii = []; jj = []; ss = [];
    for ni = 1:N
        xc = node(ni,1); yc = node(ni,2);
        phi = nodal_basis(xc,yc,h,rnode);
        mi = find(phi>0);
        mj = ni*ones(size(mi));
        ms = phi(mi);
        ii = [ii;mi]; jj = [jj;mj]; ss = [ss;ms];
    end

    rows = []; cols = []; vals = [];
    for nr = 1:Nray
        row = ii;
        col = jj + (nr-1)*N;
        phase = real(ray(jj,nr)).*rnode(ii,1) + imag(ray(jj,nr)).*rnode(ii,2);
        val = ss.*exp(1i*omega*phase);
        rows = [rows; row];
        cols = [cols; col];
        vals = [vals; val];
    end
    
    A = sparse(rows,cols,vals,rN,N*Nray);
    AT = transpose(A);
    
    b = pde.ex_u(rnode);
    u = (AT*A)\(AT*b);
    toc;
    
    fquadorder= 3;
    [err1] = Ray_FEM_solution_error(node,elem,omega,pde,ray,u,fquadorder);
    best_err(oi) = err1;
    end
    
    %% Ray-FEM
    fprintf('Ray-FEM time: \n');
    tic;
    fquadorder= 9;
    [~,~,v,~,~,err2] = Ray_FEM_IBC_1(node,elem,omega,pde,ray,fquadorder,plt);
    solu_err(oi) = err2;
    toc;
    
    err2/err1
    
end

showrate(omegas(1:test_n),solu_err)





