%% One point source outside domain problem
% 
% 1. Apply Matrix Pencil Method with sampling data on two orthogonal segments
% to estimate the low frequency ray directions
% 
% 2. Compute the high frequency equation by Ray-FEM with numerically computed
% ray directions
%

global omega xs ys a;
xs = 10; ys = 10;                 % source location
aa = 1/2;                         % computational domain [-1/2,1/2]^2

pde = MPM_data1;                  % one point source pde structure
Nray = 1;                         % one ray direction at every grid point
plt = 0;                          % plot solution or not
fquadorder = 3;                   % numerical quadrature order
solver = 'DIR';                   % direct method to solve Au=b by u=A\b  

%% record frequencies and errors
rec_omega = [25 50 100 200 400]*pi;   
rec_err1 = zeros(length(rec_omega),1);
rec_err2 = rec_err1; rec_err3 = rec_err1; rec_err4 = rec_err1; 
rec_err5 = rec_err1; rec_err6 = rec_err1; rec_err7 = rec_err1;
rec_err8 = rec_err1; rec_err9 = rec_err1; rec_err10 = rec_err1;

nn = 4;

for j = 1:nn%length(rec_omega)
    tic;
    j
    high_omega = rec_omega(j);         % high frequency
    low_omega = sqrt(high_omega);      % low frequency
    omega = high_omega;                % high frequency
    NPW = 8;                           % number of grid points per wavelength  
    h = 1/round((NPW*omega)/(2*pi));   % mesh size
    
    
    %% coarse mesh
    a = aa;
    ch = 4*h;
    [cnode,celem] = squaremesh([-a,a,-a,a],ch);
    
    
    %% solving low frequency equation on a relative larger domain
    omega = low_omega;                  % low frequency
    wl = 2*pi/omega;                    % low frequency wavelength
    rh = wl/NPW;                        % grid size on sampling segments
    
    a = aa + h*(round(wl/h) + 1);      
    [lnode,lelem] = squaremesh([-a,a,-a,a],h); 
    [lu,~,~,rel_L2_err] = Standard_FEM_IBC(lnode,lelem,omega,pde,fquadorder,solver,plt);
    eu = pde.ex_u(lnode);
    
    rec_err1(j) = norm(lu-eu,inf);      % L^{infty} error for lwo frequency solution
    rec_err2(j) = rel_L2_err;           % relative L^2 error for lwo frequency solution
    
    
    %% Ray estimation by MPM for low frequency wave
    cN = size(cnode,1);
    err3 = zeros(cN,Nray);  err4 = err3; 
    cexray = pde.ray(cnode);
    cnumray = zeros(cN,Nray);
    
    for i = 1:cN
        x0 = cnode(i,1);  y0 = cnode(i,2);
        r = -wl:rh:wl;  r = r';
        
        %% MPM sampling on the first direction
        est_ang = pi/6;
        x = r*cos(est_ang) + x0;  y = r*sin(est_ang) + y0;
        xy = [x, y];
        u = interpolation(lnode,lelem,xy,lu);
        u = u.*sqrt((x-xs).^2 + (y-ys).^2).^(1/2);
        [z] = Matrix_Pencil(u,1);
        
        dang0 = acos(1i*imag(log(z))/(1i*omega*rh));
        ang0 = [est_ang + dang0, est_ang - dang0];
        ang0 = principal_angle(ang0);
        
        
        %% MPM sampling on the second direction
        est_ang = 0;
        x = r*cos(est_ang) + x0;  y = r*sin(est_ang) + y0;
        xy = [x, y];
        u = interpolation(lnode,lelem,xy,lu);
        u = u.*sqrt((x-xs).^2 + (y-ys).^2).^(1/2);
        [z] = Matrix_Pencil(u,1);
        
        dang1 = acos(1i*imag(log(z))/(1i*omega*rh));
        ang1 = [est_ang + dang1, est_ang - dang1];
        ang1 = principal_angle(ang1);
        
        %% MPM sampling on the third direction
        est_ang =  pi/2;
        x = r*cos(est_ang) + x0;  y = r*sin(est_ang) + y0;
        xy = [x, y];
        u = interpolation(lnode,lelem,xy,lu);
        u = u.*sqrt((x-xs).^2 + (y-ys).^2).^(1/2);
        [z] = Matrix_Pencil(u, 1);
        
        dang2 = acos(1i*imag(log(z))/(1i*omega*rh));
        ang2 = [est_ang + dang2, est_ang - dang2];
        ang2 = principal_angle(ang2);
        
        %% correction
        tol = 2*pi/180;
        cnumray(i) = find_angle_3(ang0, ang1, ang2, tol);
%         cnumray(i) = find_angle(ang1, ang2);
%         err1(i) = abs(angle - cexray(i));
        err3(i) = min(abs(ang1 - cexray(i)));
        err4(i) = min(abs(ang2 - cexray(i)));     
    end
    rec_err3(j) = norm(err3(:),inf);     % L^{infty} error for MPM sampling on one direction
    rec_err4(j) = norm(err4(:),inf);     % L^{infty} error for MPM sampling on another direction
    
    
    %% interpolate the ray information on fine mesh
    a = 1/2;
    [node,elem] = squaremesh([-a,a,-a,a],h);
    N = size(node,1);
    numray = interpolation(cnode,celem,node,cnumray);
    exray = pde.ray(node);
    ray_err = exray - numray;
    rec_err5(j) = norm(ray_err(:),inf);                % L^{infty} error 
    rec_err6(j) = norm(ray_err(:))/norm(exray(:));     % relative L^2 error 

    
    %% numerical Ray-FEM
    omega = high_omega;
    numray = exp(1i*numray);
    [uh,~,~,~,rel_L2_err] = Ray_FEM_IBC_1(node,elem,omega,pde,numray,fquadorder,plt);
    ue = pde.ex_u(node);
    rec_err7(j) = norm(ue-uh,inf);        % L^{infty} error for high frequency solution
    rec_err8(j) = rel_L2_err;             % relative L^2 error for high frequency solution
    
    
    %% Exact Ray-FEM
    ray = exp(1i*exray);
    [uh,~,~,~,rel_L2_err] = Ray_FEM_IBC_1(node,elem,omega,pde,ray,fquadorder,plt);
    ue = pde.ex_u(node);
    rec_err9(j) = norm(ue-uh,inf);        % L^{infty} error for high frequency solution
    rec_err10(j) = rel_L2_err;            % relative L^2 error for high frequency solution
    toc;
end

%% print out results
fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('Point source location:  xs = %d,  ys = %d  ', xs,ys);
fprintf(['\n' '-'*ones(1,80) '\n']);

fprintf('omega/pi    ');
fprintf('  &  %.2e',rec_omega/pi');
fprintf(['\n' '-'*ones(1,80) '\n']);

fprintf('Error_1     ');  fprintf('  &  %.2e',rec_err1');  fprintf('\n\n');
fprintf('Error_2     ');  fprintf('  &  %.2e',rec_err2');  fprintf('\n\n');
fprintf('Error_3     ');  fprintf('  &  %.2e',rec_err3');  fprintf('\n\n');
fprintf('Error_4     ');  fprintf('  &  %.2e',rec_err4');  fprintf('\n\n');
fprintf('Error_5     ');  fprintf('  &  %.2e',rec_err5');  fprintf('\n\n');
fprintf('Error_6     ');  fprintf('  &  %.2e',rec_err6');  fprintf('\n\n');
fprintf('Error_7     ');  fprintf('  &  %.2e',rec_err7');  fprintf('\n\n');
fprintf('Error_8     ');  fprintf('  &  %.2e',rec_err8');  fprintf('\n\n');
fprintf('Error_9     ');  fprintf('  &  %.2e',rec_err9');  fprintf('\n\n');
fprintf('Error_10    ');  fprintf('  &  %.2e',rec_err10'); fprintf('\n\n');

fprintf(['-'*ones(1,80) '\n']);

fprintf('Error_1:  L^{infty} error of the solution for low frequency equation  \n');
fprintf('Error_2:  relative L^{2} error of the solution for low frequency equation  \n');
fprintf('Error_3:  ray angle estimation L^{infty} error with sampling direction angle 0 \n');
fprintf('Error_4:  ray angle estimation L^{infty} error with sampling direction angle pi/2 \n');
fprintf('Error_5:  ray angle estimation L^{infty} error with two orthogonal sampling direction angles 0 and pi/2 \n');
fprintf('Error_6:  ray angle estimation relative L^{2} error with two orthogonal sampling direction angles 0 and pi/2 \n');
fprintf('Error_7:  L^{infty} error of Ray-FEM solution for high frequency equation with numerical ray estimation by MPM \n');
fprintf('Error_8:  relative L^{2} error of Ray-FEM solution for high frequency equation with numerical ray estimation by MPM \n');
fprintf('Error_9:  L^{infty} error of Ray-FEM solution for high frequency equation with exact ray information \n');
fprintf('Error_10: relative L^{2} error of Ray-FEM solution for high frequency equation with exact ray information \n');

fprintf(['-'*ones(1,80) '\n']);

%% plotting
figure(11);
subplot(3,2,1);
showrate(rec_omega(1:nn),rec_err5(1:nn)');
xlabel('\omega');
ylabel('Err_{\theta, {L^{2}(\Omega)}}');


subplot(3,2,2);
showrate(rec_omega(1:nn),rec_err6(1:nn)');
xlabel('\omega');
ylabel('Err_{\theta, {L^{\infty}(\Omega)}}');

subplot(3,2,3);
showrate(rec_omega(1:nn),rec_err7(1:nn)');
xlabel('\omega');
ylabel('NR-FEM Err_{u, {L^{2}(\Omega)}}');


subplot(3,2,4);
showrate(rec_omega(1:nn),rec_err8(1:nn)');
xlabel('\omega');
ylabel('NR-FEM Err_{u, {L^{\infty}(\Omega)}}');

subplot(3,2,5);
showrate(rec_omega(1:nn),rec_err9(1:nn)');
xlabel('\omega');
ylabel('ER-FEM Err_{u, {L^{2}(\Omega)}}');


subplot(3,2,6);
showrate(rec_omega(1:nn),rec_err10(1:nn)');
xlabel('\omega');
ylabel('ER-FEM Err_{u, {L^{\infty}(\Omega)}}');

