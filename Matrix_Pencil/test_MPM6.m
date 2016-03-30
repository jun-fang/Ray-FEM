%% test 6: MPM for four point source problem with numerical data
rec_omega = [25 200 400 800 1600 3200]*pi;
rec_err1 = zeros(length(rec_omega),1);
rec_err2 = rec_err1;
global omega;
global a;
plt = 0;                   % show solution or not
fquadorder = 4; 
solver = 'DIR';            % linear system solver

           
for j = 1: 1 %length(rec_omega)
    
    omega = rec_omega(j);
    high_omega = omega;
    NPW = 80;
    h = 1/round((NPW*omega)/(2*pi));   % mesh size
    wl = 2*pi/omega;    % wavelength
    
    pde = MPM_data5;
    
    Nray = 4;
  
    % exu = @(x) (0.005*x(:,1).*x(:,2) + 0.1*x(:,2)).*exp(1i*omega*(x(:,1) + x(:,2)));
    
%     exu = @(x) exp(1i*omega*(x(:,1) + x(:,2)));
    
    a = 1/2;
    a = a + wl;
    [node,elem] = squaremesh([-a,a,-a,a],h);
    
%     omega = sqrt(omega);
    [u_std,~,~,err] = Standard_FEM_IBC(node,elem,omega,pde,fquadorder,solver,plt);

    a = 1/2;
    [cnode,celem] = squaremesh([-a,a,-a,a],8*h);
    h = 4*h;
    if omega > 50*pi
        [cnode,celem] = squaremesh([-a,a,-a,a],4*h);
    end
    if omega > 100*pi
        [cnode,celem] = squaremesh([-a,a,-a,a],8*h);
    end
    if omega > 200*pi
        [cnode,celem] = squaremesh([-a,a,-a,a],16*h);
    end 
    if omega > 400*pi
        [cnode,celem] = squaremesh([-a,a,-a,a],32*h);
    end 
    if omega > 800*pi
        [cnode,celem] = squaremesh([-a,a,-a,a],128*h);
    end 
    if omega > 1600*pi
        [cnode,celem] = squaremesh([-a,a,-a,a],256*h);
    end 
    if omega > 3200*pi
        [cnode,celem] = squaremesh([-a,a,-a,a],256*h);
    end 
    cN = size(cnode,1);
    numray = zeros(cN,Nray);
    exray = pde.ray_ang(cnode);
    
    tic;
    for i = 1:cN
        x0 = cnode(i,1);
        y0 = cnode(i,2);
        r = -wl:h:wl;
        r = r';
        est_ang = pi/3;
        x = r*cos(est_ang) + x0;
        y = r*sin(est_ang) + y0;
        xy = [x, y];
        u = interpolation(node,elem,xy,u_std);
        
        %% MPM
        [z] = Matrix_Pencil_2(u);
        [~,sz] = sort(imag(z));
        z = z(sz);  % [5/4*pi, 7/4*pi, 3/4*pi, 1/4*pi]
        
        numray(i,:) = [ est_ang - acos(log(z(4))/(1i*high_omega*h)),...
            est_ang + acos(log(z(3))/(1i*high_omega*h)),...
            est_ang + acos(log(z(1))/(1i*high_omega*h)),...
            2*pi + est_ang - acos(log(z(2))/(1i*high_omega*h)) ];
        
%         [abs(exray(i,3) - est_ang - acos(log(z(1))/(1i*omega*h)) ),...
%             abs(2*pi - (exray(i,4) - est_ang) - acos(log(z(2))/(1i*omega*h)) ),...
%             abs(exray(i,2) - est_ang - acos(log(z(3))/(1i*omega*h)) ),...
%             abs(est_ang - exray(i,1) - acos(log(z(4))/(1i*omega*h)) )];
        
        
        
        %% NMLA
%         c0 = 1;
%         Rest = 4;
%         plt = 0;
%         [angle] = NMLA_2D_2nd(x0,y0,c0,omega,Rest,node,elem,0,0,0,pde,1/8,Nray,'ex',0,plt);
%         err2(i,:) = angle - exray(i,:);
        
    end
    toc;
    
    rec_err1(j) = norm(numray(:) - exray(:))*h/(norm(exray(:))*h);    
%     rec_err2(j) = norm(err2(:))*h/(norm(exray(:))*h);
end

fprintf('%.4e   ',rec_err1');
fprintf('\n\n');
fprintf('%.4e   ',rec_err2');
fprintf('\n\n');



% figure(1);
% showrate(rec_omega, rec_err1');
% 
% figure(2);
% showrate(rec_omega, rec_err2');

