%% Test Perfect plane wave Problem: Iterative Idea
% uexact = (exp(x) + exp(y))*exp(1i*k*sqrt(2)/2*(x+y))
clear;
fileID = fopen('result2_iter.txt','w');


%% Load source data
pde = Helmholtz_data2;
fprintf(fileID,'\n\n Perfect plane wave problem: \n\n  u_ex = (exp(x) + exp(y))*exp(1i*k*sqrt(2)/2*(x+y)) \n\n');


%% Set up
plt = 0;                   % show solution or not
fquadorder = 3;            % numerical quadrature order
solver = 'HIF';            % linear system solver

rec_N = 4;                 % we test rec_N examples

% record h and omega
rec_omega = zeros(1,rec_N);
rec_h = rec_omega;

% record the error and condition number of Standard FEM
rec_S_err = rec_omega;
rec_S_cond = rec_omega;

% record the L2 error of the numerical angle estimation
rec_ang_err1 = rec_omega;   
rec_ang_err2 = rec_omega;  

% record the error and condition number of Numerical Ray-based FEM
rec_NR_err1 = rec_omega;
rec_NR_err2 = rec_omega;
rec_NR_cond1 = rec_omega;
rec_NR_cond2 = rec_omega;


% record the error and condition number of Exact Ray-based FEM
rec_ER_err = rec_omega;
rec_ER_cond = rec_omega;

% record the error and condition number of Phase-based FEM
rec_P_err = rec_omega;
rec_P_best_err = rec_omega;
rec_P_cond = rec_omega;


global omega;
global a;
lg_a = 1;
sm_a = 1/2;
high_omega = 0;
NPW = 10;
test_omega = [20*pi 35*pi 80*pi 100*pi];


for rec_i = 1: rec_N
    high_omega = test_omega(rec_i);
%     high_omega = high_omega + 20*pi;
    low_omega = sqrt(high_omega);
    h = 1/round((NPW*high_omega)/(2*pi));

    fprintf(['-'*ones(1,80) '\n']);
    fprintf(['-'*ones(1,80) '\n']);
    fprintf('\ncase %d: \nomega/(2*pi) = %d,      1/h = %d \n\n', rec_i, high_omega/(2*pi), 1/h);
    rec_omega(rec_i) = high_omega;
    rec_h(rec_i) = h;


%% Step 1: Solve the Hemholtz equation with wavenumber sqrt(k) by Standard FEM? mesh size h = 1/k
fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('Numerical Ray-based FEM: \n\n');
fprintf('Step1 \n');

omega = low_omega;
a = lg_a;

[node,elem] = squaremesh([-a,a,-a,a],h);
[u_std,~,~,~] = Standard_FEM_IBC(node,elem,omega,pde,fquadorder,solver,plt);


%% Step 2: Use NMLA to find ray directions d_c with wavenumber sqrt(k)
fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('\nStep2 \n');

N = size(node,1);
n = round(sqrt(N));
ln = n;

uu = reshape(u_std,n,n);
uux = uu;   uuy = uu;   % numerical Du

uux(:,2:n-1) = (uu(:,3:n) - uu(:,1:n-2))/(2*h);
uux(:,n) = 2*uux(:,n-1) - uux(:,n-2);
uux(:,1) = 2*uux(:,2) - uux(:,3);
ux = uux(:);


uuy(2:n-1,:) = (uu(3:n,:) - uu(1:n-2,:))/(2*h);
uuy(n,:) = 2*uuy(n-1,:) - uuy(n-2,:);
uuy(1,:) = 2*uuy(2,:) - uuy(3,:);
uy = uuy(:);


% NMLA
Rest = 4;
lnode = node;  % large domain nodes and elements
lelem = elem;

a = sm_a;
[node,elem] = squaremesh([-a,a,-a,a],h);
N = size(node,1);
n = round(sqrt(N));

ch = 1/round((NPW*low_omega)/(2*pi));
[cnode,celem] = squaremesh([-a,a,-a,a],ch);
cN = size(cnode,1);


Nray = 1;
cnumray = zeros(cN,Nray);


if (1)
    for i = 1:cN
        x0 = cnode(i,1);  y0 = cnode(i,2);
        c0 = pde.speed(cnode(i,:));
        [cnumray(i,:)] = NMLA_2D(x0,y0,c0,omega,Rest,lnode,lelem,u_std,ux,uy,pde,1/5,Nray,'num',plt);
    end
    numray = interpolation(cnode,celem,node,cnumray);
end

if (0)
    numray = zeros(N,Nray);
    for i = 1:N
        x0 = node(i,1);  y0 = node(i,2);
        c0 = pde.speed(node(i,:));
        [numray(i,:)] = NMLA_2D(x0,y0,c0,omega,Rest,node,elem,0,0,0,pde,1/5,Nray,'ex',plt);
    end
end


ray = pde.ray(node);
ray = [real(ray), imag(ray)];
ray_dir = atan2(ray(:,2),ray(:,1));
neg_index = find(ray_dir<0);
ray_dir(neg_index) = ray_dir(neg_index) + 2*pi;

diffang = numray - ray_dir;
rec_ang_err1(rec_i) = h*norm(diffang,2);

numray = exp(1i*numray);



%% Step 3: Solve the original Helmholtz equation by Ray-based FEM with ray directions d_c
fprintf(['\n' '-'*ones(1,80) '\n']);

fprintf('\nStep3 \n');

omega = high_omega;
a = sm_a;
[node,elem] = squaremesh([-a,a,-a,a],h);
% [uh,~,relerr_uh,best_err,Cond_num] = FJ_Ray_FEM_Conj(node,elem,k,pde,numray,fquadorder,opt,solver);
[uh,A,~,~,rel_L2_err] = Ray_FEM_IBC_1(node,elem,omega,pde,numray,fquadorder,plt);
rec_NR_err1(rec_i) = rel_L2_err;
rec_NR_cond1(rec_i) = condest(A);


%% Step 4: NMLA to find original ray directions d_o with wavenumber k
tic;
fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('\nStep4 \n');

a = lg_a;
ex_Du = pde.Du(lnode);
ex_u = pde.ex_u(lnode);
ux = ex_Du(:,1);
uy = ex_Du(:,2);
ux = reshape(ux,ln,ln);
uy = reshape(uy,ln,ln);
eu = reshape(ex_u,ln,ln);

uu = reshape(uh,n,n);

uux = uu;   uuy = uu;   % numerical Du
nn = round((ln-n)/2) + 1 : round((ln+n)/2); 
eu(nn,nn) = uu;

uux(:,2:n-1) = (uu(:,3:n) - uu(:,1:n-2))/(2*h);
uux(:,n) = 2*uux(:,n-1) - uux(:,n-2);
uux(:,1) = 2*uux(:,2) - uux(:,3);
ux(nn,nn) = uux;

uuy(2:n-1,:) = (uu(3:n,:) - uu(1:n-2,:))/(2*h);
uuy(n,:) = 2*uuy(n-1,:) - uuy(n-2,:);
uuy(1,:) = 2*uuy(2,:) - uuy(3,:);
uy(nn,nn) = uuy;

ux = ux(:);   uy = uy(:);  uh = eu(:);


if (1)
    for i = 1:cN
        x0 = cnode(i,1);  y0 = cnode(i,2);
        c0 = pde.speed(cnode(i,:));
        [cnumray(i,:)] = NMLA_2D(x0,y0,c0,omega,Rest,lnode,lelem,uh,ux,uy,pde,1/5,Nray,'num',plt);
    end
    numray = interpolation(cnode,celem,node,cnumray);
end

if (0)
    numray = zeros(N,Nray);
    for i = 1:N
        x0 = node(i,1);  y0 = node(i,2);
        c0 = pde.speed(node(i,:));
        [numray(i,:)] = NMLA_2D(x0,y0,c0,omega,Rest,node,elem,0,0,0,pde,1/5,Nray,'ex',plt);
    end
end

diffang = numray - ray_dir;
rec_ang_err2(rec_i) = h*norm(diffang,2);

numray = exp(1i*numray);

%% Step 5: Solve the original Helmholtz equation by Ray-based FEM with ray directions d_o
fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('\nStep5 \n');


omega = high_omega;
a = sm_a;
[node,elem] = squaremesh([-a,a,-a,a],h);
[uh,A,~,~,rel_L2_err] = Ray_FEM_IBC_1(node,elem,omega,pde,numray,fquadorder,plt);
rec_NR_err2(rec_i) = rel_L2_err;
rec_NR_cond2(rec_i) = condest(A);




%% Standard FEM
    if (1)
        fprintf('\nStandard FEM: \n');
        [~,A,~,rel_L2_err] = Standard_FEM_IBC(node,elem,omega,pde,fquadorder,solver,plt);
        rec_S_err(rec_i) = rel_L2_err;
        rec_S_cond(rec_i) = condest(A);
    end
    
    %% Exact Ray-based FEM:
    if (1) 
        fprintf('\nExact Ray-based FEM: \n');
        ray = pde.ray(node);
        [~,A,~,~,rel_L2_err] = Ray_FEM_IBC_1(node,elem,omega,pde,ray,fquadorder,plt);
        rec_ER_err(rec_i) = rel_L2_err;    
        rec_ER_cond(rec_i) = condest(A);
    end
      
    
    %% Phase-based FEM:
    if (1)
        fprintf('\nPhase-based FEM: \n');
        [~,A,~,~,rel_L2_err] = Phase_FEM_IBC(node,elem,omega,pde,fquadorder,plt);
        rec_P_err(rec_i) = rel_L2_err;
        rec_P_cond(rec_i) = condest(A);  
    end


end


%% record and print results
rec = [rec_omega;rec_h;rec_S_err;rec_S_cond;...
    rec_ang_err1;rec_ang_err2;rec_NR_err1;rec_NR_err2;rec_NR_cond1;rec_NR_cond2;...
    rec_ER_cond;rec_ER_err;rec_P_cond;rec_P_err;];

fprintf( fileID,['\n' '-'*ones(1,80) '\n']);
fprintf( fileID,'omega:                  ');
fprintf( fileID,'%1.2e  ',rec_omega );
fprintf( fileID,'\nomega/2pi:              ');
fprintf( fileID,'%1.2e  ',rec_omega/(2*pi) );
fprintf( fileID,'\n\nGrid size h:            ');
fprintf( fileID,'%1.2e  ',rec_h);
fprintf( fileID,'\n1/h:                    ');
fprintf( fileID,'%1.2e  ',1./rec_h);

fprintf( fileID,['\n' '-'*ones(1,80) '\n']);
fprintf( fileID,'\nNumerical Ray-based FEM:\n\n');
fprintf( fileID,'Angle L2 error 1:       ');
fprintf( fileID,'%1.2d  ',rec_ang_err1);
fprintf( fileID,'\n\nAngle L2 error 2:       ');
fprintf( fileID,'%1.2d  ',rec_ang_err2);
fprintf( fileID,'\n\nRelative L2 error 1:    ');
fprintf( fileID,'%1.2d  ',rec_NR_err1);
fprintf( fileID,'\n\nRelative L2 error 2:    ');
fprintf( fileID,'%1.2d  ',rec_NR_err2);
fprintf( fileID,'\n\nCondition number 1:     ');
fprintf( fileID,'%1.2d  ',rec_NR_cond1);
fprintf( fileID,'\n\nCondition number 2:     ');
fprintf( fileID,'%1.2d  ',rec_NR_cond2);


fprintf( fileID,['\n' '-'*ones(1,80) '\n']);
fprintf( fileID,'\nStandard FEM:\n\n');
fprintf( fileID,'Condition number:       ');
fprintf( fileID,'%1.2d  ',rec_S_cond);
fprintf( fileID,'\n\nRelative L2 error:      ');
fprintf( fileID,'%1.2d  ',rec_S_err);

fprintf( fileID,['\n' '-'*ones(1,80) '\n']);
fprintf( fileID,'\nExact Ray-based FEM:\n\n');
fprintf( fileID,'Condition number:       ');
fprintf( fileID,'%1.2d  ',rec_ER_cond);
fprintf( fileID,'\n\nRelative L2 error:      ');
fprintf( fileID,'%1.2d  ',rec_ER_err);

fprintf( fileID,['\n' '-'*ones(1,80) '\n']);
fprintf( fileID,'\nPhase-based FEM:\n\n');
fprintf( fileID,'Condition number:       ');
fprintf( fileID,'%1.2d  ',rec_P_cond);
fprintf( fileID,'\n\nRelative L2 error:      ');
fprintf( fileID,'%1.2d  ',rec_P_err);
fprintf( fileID,['\n' '-'*ones(1,80) '\n']);


