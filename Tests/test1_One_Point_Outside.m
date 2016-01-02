%% One Point Source Problem: outside domain
% uexact = sqrt(k)*besselh(0,1,k*sqrt((x-2)^2 + (y-2)^2))
fileID = fopen('result1.txt','w');


%% Load source data
pde = Helmholtz_data1;
fprintf(fileID,'\n\n One point source problem (outside domain): \n\n  u_ex = sqrt(k)*besselh(0,1,k*sqrt((x-2)^2 + (y-2)^2)) \n\n');


%% Set up
plt = 0;                   % show solution or not
fquadorder = 3;            % numerical quadrature order
solver = 'DIR';            % linear system solver

rec_N = 1;                 % we test rec_N examples

% record h and omega
rec_omega = zeros(1,rec_N);
rec_h = rec_omega;

% record the error and condition number of Standard FEM
rec_S_err = rec_omega;
rec_S_cond = rec_omega;

% record the error of the numerical angle estimation
rec_max_ang = rec_omega;   % maximun norm
rec_L2_ang = rec_omega;    % L2 norm

% record the error and condition number of Numerical Ray-based FEM
rec_NR_err = rec_omega;
rec_NR_cond = rec_omega;

% record the error and condition number of Exact Ray-based FEM
rec_ER_err = rec_omega;
rec_ER_cond = rec_omega;

% record the error and condition number of Phase-based FEM
rec_P_err = rec_omega;
rec_P_cond = rec_omega;


global omega;
global a;
a = 1/2;
omega = 0;
NPW = 10;

for rec_i = 1: rec_N
    
    omega = omega + 40*pi;
    h = 1/round((NPW*omega)/(2*pi));
    [node,elem] = squaremesh([-a,a,-a,a],h);

    fprintf(['-'*ones(1,80) '\n']);
    fprintf('case %d: \nomega/(2*pi) = %d,      1/h = %d \n', rec_i, omega/(2*pi), 1/h);
    rec_omega(rec_i) = omega;
    rec_h(rec_i) = h;
    
        
    %% Standard FEM
    if (0)
        fprintf('\nStandard FEM: \n');
        [~,A,~,rel_L2_err] = Standard_FEM_IBC(node,elem,omega,pde,fquadorder,solver,plt);
        rec_S_err(rec_i) = rel_L2_err;
        rec_S_cond(rec_i) = condest(A);
    end
    
    %% Numerical Ray-based FEM
    if (0)
        fprintf('\nNumerical Ray-based FEM: \n');
        Rest = 4;
        
        ch = 1/round((NPW*sqrt(omega))/(2*pi));
        [cnode,celem] = squaremesh([-a,a,-a,a],ch);
        cN = size(cnode,1);
        
        Nray = 1;
        cnumray = zeros(cN,Nray);
                
        fprintf('NMLA time: \n');
        tic;
        if (1)
            for i = 1:cN
                x0 = cnode(i,1);  y0 = cnode(i,2);
                c0 = pde.speed(cnode(i,:));
                [cnumray(i,:)] = NMLA_2D(x0,y0,c0,omega,Rest,node,elem,0,0,0,pde,1/5,Nray,'ex',plt);
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
        toc;

        ray = pde.ray(node);
        ray = [real(ray), imag(ray)];
        ray_dir = atan2(ray(:,2),ray(:,1));
        neg_index = find(ray_dir<0);
        ray_dir(neg_index) = ray_dir(neg_index) + 2*pi;
        
        diffang = numray - ray_dir;
        rec_max_ang(rec_i) = norm(diffang,inf);
        rec_L2_ang(rec_i) = h*norm(diffang,2);
        
        numray = exp(1i*numray);
        [~,A,~,~,rel_L2_err] = Ray_FEM_IBC_1(node,elem,omega,pde,numray,fquadorder,plt);
        rec_NR_err(rec_i) = rel_L2_err;
        rec_NR_cond(rec_i) = condest(A);
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
    if (0)
        fprintf('\nPhase-based FEM: \n');
        [~,A,~,~,rel_L2_err] = Phase_FEM_IBC(node,elem,omega,pde,fquadorder,plt);
        rec_P_err(rec_i) = rel_L2_err;
        rec_P_cond(rec_i) = condest(A);  
    end
      
end


%% record and print results
rec = [rec_omega;rec_h;rec_S_err;rec_S_cond;...
    rec_max_ang;rec_L2_ang;rec_NR_cond;rec_NR_err;...
    rec_ER_cond;rec_ER_err;rec_P_cond;rec_P_err;];

fprintf( fileID,['\n' '-'*ones(1,80) '\n']);
fprintf( fileID,'omega:              ');
fprintf( fileID,'%1.2e  ',rec_omega );
fprintf( fileID,'\nomega/2pi:          ');
fprintf( fileID,'%1.2e  ',rec_omega/(2*pi) );
fprintf( fileID,'\n\nGrid size h:        ');
fprintf( fileID,'%1.2e  ',rec_h);
fprintf( fileID,'\n1/h:                ');
fprintf( fileID,'%1.2e  ',1./rec_h);

fprintf( fileID,['\n' '-'*ones(1,80) '\n']);
fprintf( fileID,'\nNumerical Ray-based FEM:\n\n');
fprintf( fileID,'Maximun angle err:  ');
fprintf( fileID,'%1.2d  ',rec_max_ang);
fprintf( fileID,'\n\nAverage angle err:  ');
fprintf( fileID,'%1.2d  ',rec_L2_ang);
fprintf( fileID,'\n\nCondition number:   ');
fprintf( fileID,'%1.2d  ',rec_NR_cond);
fprintf( fileID,'\n\nRelative L2 error:  ');
fprintf( fileID,'%1.2d  ',rec_NR_err);

fprintf( fileID,['\n' '-'*ones(1,80) '\n']);
fprintf( fileID,'\nStandard FEM:\n\n');
fprintf( fileID,'Condition number:   ');
fprintf( fileID,'%1.2d  ',rec_S_cond);
fprintf( fileID,'\n\nRelative L2 error:  ');
fprintf( fileID,'%1.2d  ',rec_S_err);

fprintf( fileID,['\n' '-'*ones(1,80) '\n']);
fprintf( fileID,'\nExact Ray-based FEM:\n\n');
fprintf( fileID,'Condition number:   ');
fprintf( fileID,'%1.2d  ',rec_ER_cond);
fprintf( fileID,'\n\nRelative L2 error:  ');
fprintf( fileID,'%1.2d  ',rec_ER_err);

fprintf( fileID,['\n' '-'*ones(1,80) '\n']);
fprintf( fileID,'\nPhase-based FEM:\n\n');
fprintf( fileID,'Condition number:   ');
fprintf( fileID,'%1.2d  ',rec_P_cond);
fprintf( fileID,'\n\nRelative L2 error:  ');
fprintf( fileID,'%1.2d  ',rec_P_err);
fprintf( fileID,['\n' '-'*ones(1,80) '\n']);


%% plot
if (0)
    figure(1);
    plot(rec_omega,rec_NR_err,'r-*');
    hold on;
    plot(rec_omega,rec_ER_err,'b-o');
    hold on;
    plot(rec_omega,rec_P_err,'k-s');
    hold on;
    legend('NR-FEM','ER-FEM','P-FEM','LOCATION','Best');
    xlabel('frequency \omega');
    ylabel('Relative L2 error');
    
    
    figure(2);
    subplot(1,3,1);
    r1 = FJ_showrate(rec_omega.*rec_h.*rec_L2_ang,rec_NR_err);
    xlabel('kh*AErr');
    ylabel('Relative L2 error');
    title(['NR-FEM RErr = C_1(kh*AErr)^{' num2str(r1) '}'],'FontSize', 14);
    subplot(1,3,2);
    r2 = FJ_showrate(rec_omega.*rec_h.*rec_h,rec_ER_err);
    xlabel('kh^2');
    ylabel('Relative L2 error');
    title(['ER-FEM RErr = C_2(kh^2)^{' num2str(r2) '}'],'FontSize', 14);
    subplot(1,3,3);
    r3 = FJ_showrate(rec_omega.*rec_h.*rec_h,rec_P_err);
    xlabel('kh^2');
    ylabel('Relative L2 error');
    title(['P-FEM RErr = C_3(kh^2)^{' num2str(r3) '}'],'FontSize', 14);
    
    
    figure(3);
    r4 = FJ_showrate(rec_omega,rec_L2_ang);
    xlabel('k');
    ylabel('AErr');
    title(['NMLA AErr = C_1(k)^{' num2str(r4) '}'],'FontSize', 14);
    
end


%% showrate wrt k

if (0)
    figure(2);
    subplot(1,3,1);
    r1 = FJ_showrate(rec_omega,rec_NR_err);
    xlabel('k');
    ylabel('Relative L2 error');
    title(['NR-FEM RErr = C_1(k)^{' num2str(r1) '}'],'FontSize', 14);
    subplot(1,3,2);
    r2 = FJ_showrate(rec_omega,rec_ER_err);
    xlabel('k');
    ylabel('Relative L2 error');
    title(['ER-FEM RErr = C_2(k)^{' num2str(r2) '}'],'FontSize', 14);
    subplot(1,3,3);
    r3 = FJ_showrate(rec_omega,rec_P_err);
    xlabel('k');
    ylabel('Relative L2 error');
    title(['P-FEM RErr = C_3(k)^{' num2str(r3) '}'],'FontSize', 14);
    
end

%% showrate for ER-FEM
if (0)
    figure(4);
    r2 = FJ_showrate(rec_h,rec_ER_err);
    xlabel('h');
    ylabel('Relative L2 error');
    title(['ER-FEM RErr = C_2(h)^{' num2str(r2) '}'],'FontSize', 14);
end

%% showrate for NR-FEM AErr
if (0)
    r1 = FJ_showrate(rec_L2_ang,rec_NR_err);
    xlabel('AErr');
    ylabel('Relative L2 error');
    title(['NR-FEM RErr = C_1(AErr)^{' num2str(r1) '}'],'FontSize', 14);
end

%% showrate for NR-FEM h
if (0)
    r1 = FJ_showrate(rec_h,rec_NR_err);
    xlabel('h');
    ylabel('Relative L2 error');
    title(['NR-FEM RErr = C_1(h)^{' num2str(r1) '}'],'FontSize', 14);
end

%% showrate for NR-FEM k
if (0)
    r1 = FJ_showrate(rec_omega,rec_NR_err);
    xlabel('k');
    ylabel('Relative L2 error');
    title(['NR-FEM RErr = C_1(k)^{' num2str(r1) '}'],'FontSize', 14);
end



