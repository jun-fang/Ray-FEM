%% test 4: MPM for one point source problem, inside domain

%% What does matter for the accuracy of the angle estimation:
% 4) anlysis of the model
% 1) distance to the source point
% 2) sampling direction
% 3) measurement of the error


% rec_omega = [ 3200 ]*pi;
rec_omega = [100 200 400 800 1600]*pi;
rec_err1 = zeros(length(rec_omega),1);
rec_err2 = rec_err1;
rec_err3 = rec_err1;
global omega;
global xs;
global ys;

xs = 0;
ys = 0;

nmla = 0;

for j = 1:length(rec_omega)
    
    omega = rec_omega(j);
    NPW = 8;
    h = 1/round((NPW*omega)/(2*pi));   % mesh size
    wl = 2*pi/omega;    % wavelength
    
    pde = MPM_data4;
    Nray = 1;
    
    a = 1/2;

%     for ii = 1: 9
%         if omega > 25*pi*2^(ii-1)
%             [node,elem] = squaremesh([-a,a,-a,a],2^(ii)*2*h);
%         end
%     end
    
    [node,elem] = squaremesh([-a,a,-a,a],h);
    
    N = size(node,1);
    err1 = zeros(N,Nray);
    err2 = err1;
    err3 = err1;
    exray = pde.ray(node);
    rec_ray = 0*exray;
    count = 0;   
        
    %     tic;
    for i = 1:N
        x0 = node(i,1);
        y0 = node(i,2);
        d0 = sqrt((x0-xs)^2 + (y0-ys)^2);
        if d0 > 8*wl && d0 <= 10*wl
            r = -wl:h:wl;
            r = r';
%             est_ang =  exray(i)/1.1 + pi/12; 
            est_ang =  exray(i) + (-1)^(exray(i)>pi)*pi/30;
            rec_ray(i) = exray(i);
            x = r*cos(est_ang) + x0;
            y = r*sin(est_ang) + y0;
            xy = [x, y];
            u = pde.ex_u(xy);
            u = u.*sqrt((x-xs).^2 + (y-ys).^2).^(1/2);
            
            %% MPM
            [z] = Matrix_Pencil_2(u);
            [~,sz] = sort(imag(z));
            z = z(sz);  % [5/4*pi, 7/4*pi, 3/4*pi, 1/4*pi]
            
            
            %             if est_ang >= 0 && est_ang < pi/4
            %                 err1(i,:) = abs(abs(2*pi - exray(i) + est_ang) - acos(1i*imag(log(z(1)))/(1i*omega*h)) );
            %             end
            %
            %             if est_ang >= pi/4 && est_ang <= 2*pi
            %                 err1(i,:) = abs(abs(exray(i) - est_ang) - acos(1i*imag(log(z(1)))/(1i*omega*h)) );
            %             end
            
            err1(i,:) = abs(abs(exray(i) - est_ang) - abs(acos(1i*imag(log(z(1)))/(1i*omega*h))) );
            err2(i,:) = abs( exp(1i*omega*h*cos(exray(i) - est_ang)) - z );   
            err3(i,:) = abs(cos(exray(i) - est_ang) - 1i*imag(log(z(1)))/(1i*omega*h));
            
%             err1(i,:) = min( abs(abs(exray(i) - est_ang) - asin(sqrt(abs(1-(1i*imag(log(z(1)))/(1i*omega*h))^2))) ),...
%                 abs(abs(2*pi - exray(i) + est_ang) - asin(sqrt(1-(1i*imag(log(z(1)))/(1i*omega*h))^2)) ) );
%             
            if length(z)>1
                count = count + 1;
%                 return;
            end
            
            
            
            %% NMLA
            if (nmla)
                c0 = 1;
                Rest = d0;
                plt = 0;
                % %         uu = exu(lnode);
                % %         ux = Ux(lnode);   uy = Uy(lnode);
                % %         [angle] = NMLA_2D_2nd(x0,y0,c0,omega,Rest,lnode,lelem,uu,ux,uy,[],1/8,Nray,'num',0,plt);
                [angle] = NMLA_2D_2nd(x0,y0,c0,omega,Rest,node,elem,0,0,0,pde,1/5,Nray,'ex',0,plt);
                err2(i,:) = angle_error(angle, exray(i));
                
                [angle] = NMLA_2D_2nd(x0,y0,c0,omega,Rest,node,elem,0,0,0,pde,1/5,Nray,'ex',1,plt);
                err3(i,:) = angle_error(angle, exray(i));
            end
        end
        
    end
    %     toc;
    
%         count
    
%     rec_err1(j) = norm(err1(:))*h/(norm(exray(:))*h);
%     rec_err2(j) = norm(err2(:))*h/(norm(exray(:))*h);
%     rec_err3(j) = norm(err3(:))*h/(norm(exray(:))*h);
    
    rec_err1(j) = norm(err1(:))/norm(rec_ray);
    rec_err2(j) = norm(err1(:),inf);
    rec_err3(j) = norm(err2(:),2)*h;
end

fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('xs = %d,    ys = %d   \n\n', xs,ys);
fprintf('MPM error at ');
fprintf('Sample data angle: %d pi ', est_ang/pi);
fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('theta_1     ');
fprintf('  &  %.2e',rec_err1');
fprintf('\n\n');

fprintf('theta_2     ');
fprintf('  &  %.2e',rec_err2');
fprintf('\n\n');

fprintf('theta_3     ');
fprintf('  &  %.2e',rec_err3');
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