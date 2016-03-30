%% test 1: MPM for exact plane waves
rec_omega = [50 100 200 400 800]*pi;
rec_err1 = zeros(length(rec_omega),1);
rec_err2 = rec_err1;

for j = 1:length(rec_omega)
    
    omega = rec_omega(j);
    NPW = 6;
    h = 1/round((NPW*omega)/(2*pi));   % mesh size
    wl = 2*pi/omega;    % wavelength
    
    exu = @(x) (x(:,1).^2 + x(:,2).^2).*exp(1i*omega*(x(:,1) + x(:,2)));
    
    % exu = @(x) (0.005*x(:,1).*x(:,2) + 0.1*x(:,2)).*exp(1i*omega*(x(:,1) + x(:,2)));
    
%     exu = @(x) exp(1i*omega*(x(:,1) + x(:,2)));
    
    a = 1/2;
    
    [node,elem] = squaremesh([-a,a,-a,a],h);
    if omega > 100*pi
        [node,elem] = squaremesh([-a,a,-a,a],4*h);
    end
    if omega > 400*pi
        [node,elem] = squaremesh([-a,a,-a,a],8*h);
    end 
    N = size(node,1);
    err1 = zeros(N,1);
    err2 = err1;
    rec_z = err1;
    
    tic;
    for i = 1:N
        x0 = node(i,1);
        y0 = node(i,2);
        r = -wl:h:wl;
        r = r';
        est_ang = pi/3;
        x = r*cos(est_ang) + x0;
        y = r*sin(est_ang) + y0;
        xy = [x, y];
        u = exu(xy);
        
        [z] = Matrix_Pencil_2(u);
        rec_z(i) = z(1);
        err1(i) = abs(pi/4 + acos(log(rec_z(i))/(1i*sqrt(2)*omega*h)) - est_ang);
        err2(i) = abs(exp(1i*sqrt(2)*omega*h*cos(pi/4 - est_ang)) - rec_z(i));
    end
    toc;
    
    rec_err1(j) = norm(err1)*h/(norm(pi/4*ones(size(err1)))*h);
    
    rec_err2(j) = norm(err2)*h/(norm(rec_z)*h);
end

rec_err1
      
rec_err2

showrate(rec_omega, rec_err1');


% acos(asin(max(imag(z)))/(omega*T))