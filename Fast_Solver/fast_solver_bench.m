%% Convergence test for homogenenous case with numerical rays
function [] = fast_solver_bench(NPW, test_num)


% clear;
addpath(genpath('../../ifem/'));
addpath('../Methods/');
addpath('../NMLA/');
addpath('../Functions/')
addpath('../Plots_Prints/');
addpath(genpath('../Classes/'));
addpath('../Singularity_Treatment/');

%% Load source data
xs = 0; ys = 0;                     % point source location
epsilon = 1/(2*pi);                % cut-off parameter   
speed = @(x) ones(size(x,1),1);      % medium speed


%% Set up
plt = 0;                   % show solution or not
fquadorder = 3;            % numerical quadrature order
Nray = 1;                  % one ray direction at each grid node
sec_opt = 0;               % NMLA second order correction or not

% NPW = 4;                   % number of points per wavelength
% test_num = 6;              % we test test_num examples

% frequency
high_omega = [120 160 200 240 320 480 640 800 960]*pi;
low_omega = 2*sqrt(high_omega); 

% error
low_max_rayerr = 0*high_omega;     % L_inf ray error of low-freq waves
low_l2_rayerr = 0*high_omega;      % L_2 ray error of low-freq waves
high_max_rayerr = 0*high_omega;    % L_inf ray error of high-freq waves
high_l2_rayerr = 0*high_omega;     % L_2 ray error of high-freq waves

max_err = 0*high_omega;            % L_inf error of the numerical solution
rel_max_err = 0*high_omega;        % relative L_inf error of the numerical solution
l2_err = 0*high_omega;             % L_2 error of the numerical solution
rel_l2_err = 0*high_omega;         % relative L_2 error of the numerical solution

% wavelength
high_wl = 2*pi./high_omega;
low_wl = 2*pi./low_omega;

% mesh size
fh = 1./(NPW*round(high_omega/(2*pi)));      % fine mesh size
ch = 1./(10*round(sqrt(2./fh)/10));

% ch = 2*fh;
% ch = 1./(10*round(low_omega/(2*pi)));        % coarse mesh size
% ch = 1./(NPW*round(1./sqrt(fh)/10)*10);
% ch = fh.*ceil(ch./fh);

% width of PML
high_wpml = 0.07*ones(size(high_omega));
low_wpml = 0.25*ones(size(high_omega));

npml_high = round(high_wpml./fh);
npml_low  = round(low_wpml./ch);

high_wpml = npml_high.*fh;
low_wpml = npml_low.*ch;

% number of pml points at the interfaces 
npml_interior = 10*ones(size(high_omega));


% high_wpml = 8*high_wl(1)*ones(size(high_omega)); %fh.*ceil(high_wl./fh);
% low_wpml = ch.*ceil(low_wl(1)./ch);

% We need to define different PML's for the different cases 


%% Generate the domain sizes
sd = 1/2;
Rest = 0.4654; 2*epsilon;           % estimate of the distance to the source point

high_r = NMLA_radius(high_omega(1),Rest);
md = sd + high_r + high_wpml;
md = ceil(md*10)/10;      % middle domain size 

% Rest = sqrt(2)*md;
low_r = NMLA_radius(low_omega(1),Rest);
ld = md + low_r + low_wpml;
ld = ceil(ld*10)/10;      % large domain size



%% Tests
tstart = tic;
for ti = 1: test_num
    omega = high_omega(ti);
    h = fh(ti);  h_c = ch(ti);
    fprintf(['-'*ones(1,80) '\n']);
    fprintf('\nCase %d: \nomega/(2*pi) = %d,   1/h = %d   1/h_c = %d,  NPW = %d ',...
        ti, round(omega/(2*pi)), 1/h, 1/h_c, NPW);
    
    
    %% Step 1: Solve the Hemholtz equation with the same source but with a relative low frequency sqrt(\omega) by Standard FEM, mesh size \omega*h = constant
    fprintf(['\n' '-'*ones(1,80) '\n']);
    fprintf('Step1: S-FEM, low frequency \n');
    tic;
    omega = low_omega(ti);              % low frequency
    a = ld(ti);                         % large domain 
    wpml = low_wpml(ti);                % width of PML
    npml = npml_low(ti);
    sigmaMax = 25/wpml;                 % Maximun absorbtion
    [lnode,lelem] = squaremesh([-a,a,-a,a],h);
    nSub = round(omega / (10*pi));      % number of subdomains
    
    % smooth part
    %A = assemble_Helmholtz_matrix_SFEM(lnode,lelem,omega,wpml,sigmaMax,speed,fquadorder);
    b = assemble_RHS_SFEM_with_ST(lnode,lelem,xs,ys,omega,wpml,sigmaMax,epsilon,fquadorder);

    % uing the polarized traces solver 
    
    % setting the width of the pml layers 
    npmlInt = npml_interior(ti);
    wpmlInt = npmlInt*h;
    
    M0 = model;
    M0.init(lnode,lelem,omega,wpml,sigmaMax,speed,fquadorder);

    % extracting the interior degrees of freedom
    fInt = b(M0.freeNode);
    

    
    
    %%
    x = unique(M0.node(:,1));
    y = unique(M0.node(:,2));

    
    % gotta be careful with the PML (they have to be outside the
    % interior domain
    xInt = x(npml:end-npml);
    nInt = length(xInt);
    
    % number of points per subdomains
    nSubInt = round(nInt/nSub);
    
    % the limits for the domaind decomposition
    nIndLim = round(linspace(1,nInt-1,nSub+1));
    
    % indices at which each subdomain starts (and finishes)
    indn = npml + nIndLim(2:end);
    ind1 = npml + nIndLim(1:end-1)+1;
    ind1(1) = ind1(1) -1;
    
    % last and first coordinate
    xn = x(indn);
    x1 = x(ind1);
    
    nodeArray = {};
    elemArray = {};
    for ii = 1:nSub
        if ii == 1 % first subarray 
            [nodeArray{ii},elemArray{ii}] = squaremesh([x(ind1(ii)-npml),...
                x(indn(ii)+npmlInt),...
                -a,a],h);
        elseif ii == nSub
            [nodeArray{ii},elemArray{ii}] = squaremesh([x(ind1(ii)-npmlInt),...
                x(indn(ii)+npml),...
                -a,a],h);
        else
            [nodeArray{ii},elemArray{ii}] = squaremesh([x(ind1(ii)-npmlInt),...
                x(indn(ii)+npmlInt),...
                -a,a],h);
        end
    end
    
    % we define the physical location of the interface boundaries
    separatorN = xn(1:end-1);
    separator1 = x1(2:end);
    
    % building the array of objects
    MArray = {};
    SubArray = {};
        
    for ii = 1:nSub
        % we will need to make sure that the wpml has the correct width
        MArray{ii} = model;
        % we need to be carefull when defining the PML
        if ii == 1
            wpmlvec = [wpml, wpmlInt-2*h, wpml, wpml];
        elseif ii == nSub
            wpmlvec = [wpmlInt-2*h, wpml, wpml, wpml];
        else
            wpmlvec = [wpmlInt-2*h, wpmlInt-2*h, wpml, wpml];
        end
        MArray{ii}.init(nodeArray{ii},elemArray{ii},omega,...
            wpmlvec,h/NPW,speed,fquadorder);
        
    end
    
    %% we need to define all the local (and global) local Indices
    % We start by defining the indices associated to each boundary (locally)
    indxn = {};
    indxnp = {};
    indx0  = {};
    indx1  = {};
    
    for ii = 1:nSub
        if ii ~= 1
            indx0{ii} = find( abs(MArray{ii}.node(MArray{ii}.freeNode,1) - separatorN(ii-1)) < h/10 );
            indx1{ii} = find( abs(MArray{ii}.node(MArray{ii}.freeNode,1) - separator1(ii-1)) < h/10 );
        end
        
        if ii ~= nSub
            indxn{ii} = find( abs(MArray{ii}.node(MArray{ii}.freeNode,1) - separatorN(ii)) < h/10 );
            indxnp{ii}= find( abs(MArray{ii}.node(MArray{ii}.freeNode,1) - separator1(ii)) < h/10 );
        end
    end
    indxn{nSub} = [];
    indxnp{nSub} = [];
    
    % building the the local-global indices
    
    indIntGlobal = {};
    for ii = 1:nSub
        if ii == 1
            indIntGlobal{ii} = find(M0.node(M0.freeNode,1) <= separatorN(ii) );
        elseif ii == nSub
            indIntGlobal{ii} = find(M0.node(M0.freeNode,1) > separator1(ii-1) );
        else
            indIntGlobal{ii} = find((M0.node(M0.freeNode,1) <= separatorN(ii) + h/10 ).*  ...
                (M0.node(M0.freeNode,1) > separator1(ii-1) - h/10 ));
        end
    end
    
    
    % building the local-local indices
    
    indIntLocal = {};
    for ii = 1:nSub
        if ii == 1
            indIntLocal{ii} = find(MArray{ii}.node(MArray{ii}.freeNode,1) <= separatorN(ii) );
        elseif ii == nSub
            indIntLocal{ii} = find(MArray{ii}.node(MArray{ii}.freeNode,1) > separator1(ii-1) );
        else
            indIntLocal{ii} = find((MArray{ii}.node(MArray{ii}.freeNode,1) <= separatorN(ii)+h/10 ).*  ...
                (MArray{ii}.node(MArray{ii}.freeNode,1) > separator1(ii-1)-h/10 ));
        end
    end
    
    
    for ii = 1:nSub
        if size(indIntLocal{ii},1) ~= size(indIntGlobal{ii},1)
            fprintf('Mismatch of number of indices ot %d th layer \n', ii);
        end
    end
    
    
    
    %% building the subdomains from the models and the local indices
    for ii = 1:nSub
        if ii == 1
            position = 'N';
        elseif ii == nSub
            position = 'S';
        else
            position = 'C'; %center
        end
        
        SubArray{ii} = subdomain;
        SubArray{ii}.init_Model(MArray{ii},indIntLocal{ii}, indIntGlobal{ii}, indxn{ii}, indxnp{ii}, indx0{ii}, indx1{ii}, position)
    end
    
    %% Factorizing the local matrices
    tic();
    for ii = 1:nSub
        SubArray{ii}.LU_factorization();
    end
    toc();
    
    %% bulding the global array
    DDM = layer_stripping;
    
    DDM.init_model(SubArray,M0, nSub);
    
    %% Solving local problems
    
    uArray =  DDM.solve_local(fInt);
    
    fGamma = DDM.formRHS(uArray);
    
    %% Testing the application of the matrix M
    %
    % DDM.MBase.LU_factorization();
    % uGlobal = DDM.MBase.solveInt(fInt);
    % uGamma = DDM.extractGlobalTraces(uGlobal);
    %
    % res =  DDM.applyM(uGamma) + fGamma;
    %
    % norm(res)
    
    %% testing the application of the  matrix MM
    fGammaPol = DDM.formPolarizedRHS(uArray);
    
    MMu = DDM.applyMpolarized(fGammaPol );
    
    DDM.initPermutationMatrix();
    MM = @(u) DDM.applyMpolarized(u);
    
    %% testing the Dinv
    LL = @(u) DDM.applyGSinv(u);
    
    tic();
    uGmres = gmres(MM, -fGammaPol,[], 1e-9, 300, LL );
    toc()
    %% testing that the solution coincides
    
    u = uGmres(1:end/2) + uGmres(end/2+1:end);
    
    res = DDM.applyMpolarized(uGmres) + fGammaPol;
    norm(res)
    
    
    %% reconstructing the solution within the subdomain
    UsolArray = DDM.solveLocal(fInt, u);
    U = DDM.reconstruct(UsolArray);
    
    u_std = zeros(size(f,1),1);
    u_std(M0.freeNode) = U;
    
    % plotting the global solution
    if plt == 1
        M0.showresult(real(u_std));
    end
    %% building the full solution from the smooth and singular part
    
    x = lnode(:,1); y = lnode(:,2);
    rr = sqrt((x-xs).^2 + (y-ys).^2);
    ub = 1i/4*besselh(0,1,omega*rr);
    cf = cutoff(epsilon,2*epsilon,nodeLow,xs,ys);
    
    % low frequency solution: smooth + singularity
    u_low = u_std + ub.*cf;
    toc;
    
    if plt == 1
        M0.showresult(real(u_low));
    end
    

    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [~,~,isBdNode] = findboundary(lelem);
    freeNode = find(~isBdNode);
    lN = size(lnode,1);        u_std = zeros(lN,1);
    u_std(freeNode) = A(freeNode,freeNode)\b(freeNode);
    
    % singular part
    x = lnode(:,1); y = lnode(:,2);
    rr = sqrt((x-xs).^2 + (y-ys).^2);
    ub = 1i/4*besselh(0,1,omega*rr);
    cf = cutoff(epsilon,2*epsilon,lnode,xs,ys);
    
    % low frequency solution: smooth + singularity
    u_low = u_std + ub.*cf;
    toc;
    
    
    %% Step 2: Use NMLA to find ray directions d_c with low frequency sqrt(\omega)
    fprintf([ '-'*ones(1,80) '\n']);
    fprintf('Step2: NMLA, low frequency \n');
    
    % compute numerical derivatives 
    [ux,uy] = num_derivative(u_low,h,2);
    
    a = md(ti);
    [mnode,melem] = squaremesh([-a,a,-a,a],h);
    
    [cnode,celem] = squaremesh([-a,a,-a,a],h_c);
    cN = size(cnode,1);
    cnumray_angle = zeros(cN,Nray);
    
    % NMLA
    tic;
% mesh size
fh = 1./(NPW*round(high_omega/(2*pi)));      % fine mesh size
ch = 1./(10*round(sqrt(2./fh)/10));

% ch = 2*fh;
% ch = 1./(10*round(low_omega/(2*pi)));        % coarse mesh size
% ch = 1./(NPW*round(1./sqrt(fh)/10)*10);
% ch = fh.*ceil(ch./fh);

% width of PML
high_wpml = 0.07*ones(size(high_omega));
low_wpml = 0.19*ones(size(high_omega));
% high_wpml = 8*high_wl(1)*ones(size(high_omega)); %fh.*ceil(high_wl./fh);
% low_wpml = ch.*ceil(low_wl(1)./ch);

    for i = 1:cN
        x0 = cnode(i,1);  y0 = cnode(i,2);
        r0 = sqrt((x0-xs)^2 + (y0-ys)^2);
        c0 = speed(cnode(i,:));
        if r0 > (2*epsilon - 3*h_c)
%             Rest = r0;
            [cnumray_angle(i)] = NMLA(x0,y0,c0,omega,Rest,lnode,lelem,u_low,ux,uy,[],1/5,Nray,'num',sec_opt,plt);
        else
            cnumray_angle(i) = ex_ray([x0,y0],xs,ys,0);
        end
    end
    cnumray = exp(1i*cnumray_angle);
    numray1 = interpolation(cnode,celem,mnode,cnumray);
    numray1 = numray1./abs(numray1);
    toc;
    
    % compute the ray errors
    exray = ex_ray(mnode,xs,ys,1);
    mr = sqrt((mnode(:,1)-xs).^2 + (mnode(:,2)-ys).^2);
    numray1 = numray1.*(mr > 2*epsilon) + exray.*(mr <= 2*epsilon);
    rayerr1 = numray1 - exray;
    low_max_rayerr(ti) = norm(rayerr1,inf);
    low_l2_rayerr(ti) = norm(rayerr1)*h/(norm(exray.*(mr > 2*epsilon))*h);
    
    
    %% Step 3: Solve the original Helmholtz equation by Ray-based FEM with ray directions d_c
    fprintf(['-'*ones(1,80) '\n']);
    fprintf('Step3: Ray-FEM, high frequency \n');
   
    tic;
    omega = high_omega(ti);
    wpml = high_wpml(ti);                % width of PML
    sigmaMax = 25/wpml;                 % Maximun absorbtion
    ray = numray1;

    % smooth part
    option = 'homogeneous';
    A = assemble_Helmholtz_matrix_RayFEM(mnode,melem,omega,wpml,sigmaMax,speed,ray,fquadorder);
    b = assemble_RHS_RayFEM_with_ST(mnode,melem,xs,ys,omega,epsilon,wpml,sigmaMax,ray,speed,fquadorder,option);
    uh = RayFEM_direct_solver(mnode,melem,A,b,omega,ray,speed);
    
    % singularity part
    x = mnode(:,1); y = mnode(:,2);
    rr = sqrt((x-xs).^2 + (y-ys).^2);
    ub = 1i/4*besselh(0,1,omega*rr);
    cf = cutoff(epsilon,2*epsilon,mnode,xs,ys);
    
    % smooth + singularity
    uh1 = uh + ub.*cf;
    toc;
    
    
    
    %% Step 4: NMLA to find original ray directions d_o with wavenumber k
    fprintf(['-'*ones(1,80) '\n']);
    fprintf('Step4: NMLA, high frequency \n');
    
    % compute numerical derivatives
    [ux,uy] = num_derivative(uh1,h,2);
    
    a = sd;
    [node,elem] = squaremesh([-a,a,-a,a],h);
    
    [cnode,celem] = squaremesh([-a,a,-a,a],h_c);
    cN = size(cnode,1);
    cnumray_angle = zeros(cN,Nray);
    
    
    % NMLA
    tic;
    for i = 1:cN
        x0 = cnode(i,1);  y0 = cnode(i,2);
        r0 = sqrt((x0-xs)^2 + (y0-ys)^2);
        c0 = speed(cnode(i,:));
        if r0 > (2*epsilon - 3*h_c)
%             Rest = r0;
            [cnumray_angle(i)] = NMLA(x0,y0,c0,omega,Rest,mnode,melem,uh1,ux,uy,[],1/5,Nray,'num',sec_opt,plt);
        else
            cnumray_angle(i) = ex_ray([x0,y0],xs,ys,0);
        end
    end
    cnumray = exp(1i*cnumray_angle);
    numray2 = interpolation(cnode,celem,node,cnumray);
    numray2 = numray2./abs(numray2);
    toc;
    
    % compute the ray errors
    exray = ex_ray(node,xs,ys,1);
    sr = sqrt((node(:,1)-xs).^2 + (node(:,2)-ys).^2);
    numray2 = numray2.*(sr > 2*epsilon) + exray.*(sr <= 2*epsilon);
    rayerr2 = numray2 - exray;
    high_max_rayerr(ti) = norm(rayerr2,inf);
    high_l2_rayerr(ti) = norm(rayerr2)*h/(norm(exray.*(sr > 2*epsilon))*h);    
    
    
    %% Step 5: Solve the original Helmholtz equation by Ray-based FEM with ray directions d_o
    fprintf(['-'*ones(1,80) '\n']);
    fprintf('Step5: Ray-FEM, high frequency \n');
    tic;
    
    % Ray-FEM solution
    omega = high_omega(ti);
    wpml = high_wpml(ti);                % width of PML
    sigmaMax = 25/wpml;                 % Maximun absorbtion
    ray = numray2;
    
    option = 'homogeneous';
    A = assemble_Helmholtz_matrix_RayFEM(node,elem,omega,wpml,sigmaMax,speed,ray,fquadorder);
    b = assemble_RHS_RayFEM_with_ST(node,elem,xs,ys,omega,epsilon,wpml,sigmaMax,ray,speed,fquadorder,option);
    u = RayFEM_direct_solver(node,elem,A,b,omega,ray,speed);
    toc;
    
    % Excat solution 
    x = node(:,1); y = node(:,2);
    rr = sqrt((x-xs).^2 + (y-ys).^2);
    
    ub = 1i/4*besselh(0,1,omega*rr);
    cf = cutoff(epsilon,2*epsilon,node,xs,ys);
    uex = (1-cf).*ub;
    uex(rr <= epsilon) = 0;
   
    % Errors
    du = u - uex;
    idx = find( ~( (x <= max(x)-wpml).*(x >= min(x)+wpml)...
        .*(y <= max(y)-wpml).*(y >= min(y)+wpml) ) ); % index on PML
    du(idx) = 0;  uex(idx) = 0;
    
    max_err(ti) = norm(du,inf);
    rel_max_err(ti) = norm(du,inf)/norm(uex,inf);
    l2_err(ti) = norm(du)*h;
    rel_l2_err(ti) = norm(du)/norm(uex);
       
end

totaltime = toc(tstart);
fprintf('\n\nTotal running time: % d minutes \n', totaltime/60);



%% save output 
nameFile = strcat('resutls_2_HomNumRay_NPW_', num2str(NPW), '.mat');
save(nameFile, 'rel_l2_err', 'NPW', 'high_omega', 'test_num');



%% plots
% figure(1);
% subplot(2,2,1);
% show_convergence_rate(high_omega(1:test_num),low_max_rayerr(1:test_num),'omega','low max');
% subplot(2,2,2);
% show_convergence_rate(high_omega(1:test_num),low_l2_rayerr(1:test_num),'omega','low l2');
% subplot(2,2,3);
% show_convergence_rate(high_omega(1:test_num),high_max_rayerr(1:test_num),'omega','high max');
% subplot(2,2,4);
% show_convergence_rate(high_omega(1:test_num),high_l2_rayerr(1:test_num),'omega','high l2');
% 
% figure(2);
% subplot(2,2,1);
% show_convergence_rate(high_omega(1:test_num),max_err(1:test_num),'omega','max err');
% subplot(2,2,2);
% show_convergence_rate(high_omega(1:test_num),l2_err(1:test_num),'omega','L2 err');
% subplot(2,2,3);
% show_convergence_rate(high_omega(1:test_num),rel_max_err(1:test_num),'omega','Rel max ');
% subplot(2,2,4);
% show_convergence_rate(high_omega(1:test_num),rel_l2_err(1:test_num),'omega','Rel L2 ');


%% print results
fprintf( ['\n' '-'*ones(1,80) '\n']);
fprintf( 'omega:                   ');
fprintf( '&  %.2e  ',high_omega );
fprintf( '\nomega/2pi:               ');
fprintf( '&  %.2e  ',high_omega/(2*pi) );
fprintf( '\n\nGrid size h:             ');
fprintf( '&  %.2e  ',fh);
fprintf( '\n1/h:                     ');
fprintf( '&  %.2e  ',1./fh);

fprintf( ['\n' '-'*ones(1,80) '\n']);
fprintf( 'Low max ray error:       ');
fprintf( '&  %1.2d  ',low_max_rayerr);
fprintf( '\n\nLow rel l2 ray error:    ');
fprintf( '&  %1.2d  ',low_l2_rayerr);
fprintf( '\n\nHigh max ray error:      ');
fprintf( '&  %1.2d  ',high_max_rayerr);
fprintf( '\n\nHigh rel l2 ray error:   ');
fprintf( '&  %1.2d  ',high_l2_rayerr);
fprintf( '\n\nMax error:               ');
fprintf( '&  %1.2d  ',max_err);
fprintf( '\n\nRelative max error:      ');
fprintf( '&  %1.2d  ',rel_max_err);
fprintf( '\n\nL2 error:                ');
fprintf( '&  %1.2d  ',l2_err);
fprintf( '\n\nRelative L2 error:       ');
fprintf( '&  %1.2d  ',rel_l2_err);


fprintf( ['\n' '-'*ones(1,80) '\n']);


