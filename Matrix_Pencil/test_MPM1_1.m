%% test 1: MPM for one perfect plane wave
% u(x,y) = ( x+y ) * exp(i*omega*sqrt(2)/2*(x+y))
% poles z1 = exp(h*cos(ang_est) + 1i*omega*h*cos(pi/4 - ang_est))
%       z2 = exp(h*sin(ang_est) + 1i*omega*h*cos(pi/4 - ang_est)) 


rec_omega = [100 200 400 800 1600]*pi;
rec_err1 = zeros(length(rec_omega),1);
rec_err2 = rec_err1;     
rec_err3 = rec_err1;           
        
global omega;

for j = 1:length(rec_omega)
    
    omega = rec_omega(j);
    NPW = 6;
    h = 1/round((NPW*omega)/(2*pi));   % mesh size
    wl = 2*pi/omega;    % wavelength
    
    exu = @(x) (x(:,1) + x(:,2)).*exp(1i*omega*sqrt(2)/2*(x(:,1) - x(:,2)));    
%     pde = Helmholtz_data2;

    a = 1/2;
    
    [node,elem] = squaremesh([-a,a,-a,a],h);
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
        [node,elem] = squaremesh([-a,a,-a,a],256*h);
    end 
    N = size(node,1);
    err1 = zeros(N,1);
    err2 = err1;
    err3 = err1;
    exray = pi/4*ones(N,1);
    
    tic;
    for i = 1:N
        x0 = node(i,1);
        y0 = node(i,2);
        r = -wl:h:wl;
        r = r';
        est_ang = pi/4;
        x = r*cos(est_ang) + x0;
        y = r*sin(est_ang) + y0;
        xy = [x, y];
        u = exu(xy);
        count = 0;
        
        %% Matrix Pencil
        [z] = Matrix_Pencil_2(u);
        [~,sz] = sort(imag(z));
        z = z(sz);
        if length(z) > 1
            count = count + 1;
        end
        rec_z = z(1);
        err1(i) = abs(pi/4 - acos(1i*imag(log(rec_z))/(1i*omega*h)) + est_ang);
        
    end
    toc;
    count 
    rec_err1(j) = norm(err1)*h/(norm(exray(:))*h);
    
    rec_err2(j) = norm(err2)*h/(norm(exray(:))*h);
    
    rec_err3(j) = norm(err3)*h/(norm(exray(:))*h);

end

fprintf('%.4e   ',rec_err1');
fprintf('\n\n');
fprintf('%.4e   ',rec_err2');
fprintf('\n\n');
fprintf('%.4e   ',rec_err3');
fprintf('\n\n');

% showrate(rec_omega, rec_err1');


