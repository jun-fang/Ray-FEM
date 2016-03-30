%% test 4: MPM for one point source problem
rec_omega = [100 200 400 800 1600 3200]*pi;
rec_err1 = zeros(length(rec_omega),1);
rec_err2 = rec_err1;
rec_err3 = rec_err1;
global omega;
global xs;
global ys;

xs = 1;
ys = 1;

nmla = 0;

for j = 1: length(rec_omega)
    
    omega = rec_omega(j);
    NPW = 6;
    h = 1/round((NPW*omega)/(2*pi));   % mesh size
    wl = 2*pi/omega;    % wavelength
    
%     exu = @(x) ( (exp(x(:,1)) + exp(x(:,2))).*exp(1i*omega*sqrt(2)/2*(x(:,1) + x(:,2))) ...
%         - 10*(x(:,1) + x(:,2)).*exp(1i*omega*sqrt(2)/2*(x(:,1) - x(:,2))) );
%     
%     Ux = @(x) ( (exp(x(:,1)) + 1i*omega*sqrt(2)/2*(exp(x(:,1)) + exp(x(:,2))))...
%         .*exp(1i*omega*sqrt(2)/2*(x(:,1) + x(:,2))) ...
%         - 10*(1 + 1i*omega*sqrt(2)/2*(x(:,1) + x(:,2)))...
%         .*exp(1i*omega*sqrt(2)/2*(x(:,1) - x(:,2))) );
%     
%     Uy = @(x) ( (exp(x(:,2)) + 1i*omega*sqrt(2)/2*(exp(x(:,1)) + exp(x(:,2))))...
%         .*exp(1i*omega*sqrt(2)/2*(x(:,1) + x(:,2))) ...
%         - 10*(1 - 1i*omega*sqrt(2)/2*(x(:,1) + x(:,2)))...
%         .*exp(1i*omega*sqrt(2)/2*(x(:,1) - x(:,2))) );
    
    pde = MPM_data4;
    
    Nray = 1;
  
    % exu = @(x) (0.005*x(:,1).*x(:,2) + 0.1*x(:,2)).*exp(1i*omega*(x(:,1) + x(:,2)));
    
%     exu = @(x) exp(1i*omega*(x(:,1) + x(:,2)));
    
    a = 1/2;
    b = 2*a;
    [node,elem] = squaremesh([-a,a,-a,a],h);
%     [lnode,lelem] = squaremesh([-b,b,-b,b],h);
    if omega > 50*pi
        [node,elem] = squaremesh([-a,a,-a,a],4*h);
    end
    if omega > 100*pi
        [node,elem] = squaremesh([-a,a,-a,a],8*h);
    end
    if omega > 200*pi
        [node,elem] = squaremesh([-a,a,-a,a],16*h);
    end 
    if omega > 400*pi
        [node,elem] = squaremesh([-a,a,-a,a],32*h);
    end 
    if omega > 800*pi
        [node,elem] = squaremesh([-a,a,-a,a],128*h);
    end 
    if omega > 1600*pi
        [node,elem] = squaremesh([-a,a,-a,a],256*h);
    end 
    if omega > 3200*pi
        [node,elem] = squaremesh([-a,a,-a,a],512*h);
    end 
    N = size(node,1);
    err1 = zeros(N,Nray);
    err2 = err1;
    err3 = err1;
    exray = pde.ray(node);
    count = 0;
    
%     tic;
    for i = 1:N
        x0 = node(i,1);
        y0 = node(i,2);
        r = -wl:h:wl;
        r = r';
        est_ang = 2.5*pi/4;
        x = r*cos(est_ang) + x0;
        y = r*sin(est_ang) + y0;
        xy = [x, y];
        u = pde.ex_u(xy);
        
        %% MPM
        [z] = Matrix_Pencil_2(u);
        [~,sz] = sort(imag(z));
        z = z(sz);  % [5/4*pi, 7/4*pi, 3/4*pi, 1/4*pi]
        if length(z)>1
            count = count +1;
        end
        
        if est_ang >= 0 && est_ang < pi/4
            err1(i,:) = abs(abs(2*pi - exray(i) + est_ang) - acos(1i*imag(log(z(1)))/(1i*omega*h)) );
        end
        
        if est_ang >= pi/4 && est_ang <= 2*pi
            err1(i,:) = abs(abs(exray(i) - est_ang) - acos(1i*imag(log(z(1)))/(1i*omega*h)) );
        end
%         if est_ang > 5*pi/4 && est_ang <= 2*pi
%             err1(i,:) = abs(abs(exray(i) - est_ang) - acos(1i*imag(log(z(1)))/(1i*omega*h)) );
%         end

        
        
        
        %% NMLA
        if (nmla)
            c0 = 1;
            Rest = 4;
            plt = 0;
            % %         uu = exu(lnode);
            % %         ux = Ux(lnode);   uy = Uy(lnode);
            % %         [angle] = NMLA_2D_2nd(x0,y0,c0,omega,Rest,lnode,lelem,uu,ux,uy,[],1/8,Nray,'num',0,plt);
            [angle] = NMLA_2D_2nd(x0,y0,c0,omega,Rest,node,elem,0,0,0,pde,1/5,Nray,'ex',0,plt);
            err2(i,:) = angle - exray(i);
            
            [angle] = NMLA_2D_2nd(x0,y0,c0,omega,Rest,node,elem,0,0,0,pde,1/5,Nray,'ex',1,plt);
            err3(i,:) = angle - exray(i);
        end

    end
%     toc;

%     count

    rec_err1(j) = norm(err1(:))*h/(norm(exray(:))*h);    
    rec_err2(j) = norm(err2(:))*h/(norm(exray(:))*h);
    rec_err3(j) = norm(err3(:))*h/(norm(exray(:))*h);
end

fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('xs = %d,    ys = %d   \n\n', xs,ys);
fprintf('MPM error at ');
fprintf('Sample data angle: %d pi ', est_ang/pi);
fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('theta_1     ');
fprintf('  &  %.2e',rec_err1');
fprintf('\n\n');

if (nmla)
    fprintf('NMLA error:\n');
    fprintf('theta_1    ');
    fprintf('  &  %.2e',rec_err2');
    fprintf('\n\n');
    
    fprintf('corrected NMLA error:\n');
    fprintf('theta_1 pi/4    ');
    fprintf('  &  %.2e',rec_err3');
    fprintf('\n\n');
end


% figure(1);
% showrate(rec_omega, rec_err1');
% 
% figure(2);
% showrate(rec_omega, rec_err2');
% 
% figure(3);
% showrate(rec_omega, rec_err3');