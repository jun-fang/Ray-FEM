%% test 2: MPM for four perfect plane waves
rec_omega = [100 200 400 800 1600]*pi;
rec_err1 = zeros(length(rec_omega),1);   
rec_err2 = rec_err1;         

for j = 1: length(rec_omega)
    
    omega = rec_omega(j);
    NPW = 6;
    h = 1/round((NPW*omega)/(2*pi));   % mesh size
    wl = 2*pi/omega;    % wavelength
    
    exu = @(x) ( 3*exp(1i*omega*sqrt(2)/2*(x(:,1) + x(:,2))) ...
        - 2*exp(1i*omega*sqrt(2)/2*(x(:,1) - x(:,2)))...
        + 4*exp(1i*omega*sqrt(2)/2*(-x(:,1) - x(:,2)))...
        + 0.8*exp(1i*omega*sqrt(2)/2*(-x(:,1) + x(:,2))) );
    
    Ux = @(x) ( 3*1i*omega*sqrt(2)/2*exp(1i*omega*sqrt(2)/2*(x(:,1) + x(:,2))) ...
        - 2*1i*omega*sqrt(2)/2*exp(1i*omega*sqrt(2)/2*(x(:,1) - x(:,2)))...
        - 4*1i*omega*sqrt(2)/2*exp(1i*omega*sqrt(2)/2*(-x(:,1) - x(:,2)))...
        - 0.8*1i*omega*sqrt(2)/2*exp(1i*omega*sqrt(2)/2*(-x(:,1) + x(:,2))) );
    
    Uy = @(x) ( 3*1i*omega*sqrt(2)/2*exp(1i*omega*sqrt(2)/2*(x(:,1) + x(:,2))) ...
        + 2*1i*omega*sqrt(2)/2*exp(1i*omega*sqrt(2)/2*(x(:,1) - x(:,2)))...
        - 4*1i*omega*sqrt(2)/2*exp(1i*omega*sqrt(2)/2*(-x(:,1) - x(:,2)))...
        + 0.8*1i*omega*sqrt(2)/2*exp(1i*omega*sqrt(2)/2*(-x(:,1) + x(:,2))) );
    
    Nray = 4;
  
    % exu = @(x) (0.005*x(:,1).*x(:,2) + 0.1*x(:,2)).*exp(1i*omega*(x(:,1) + x(:,2)));
    
%     exu = @(x) exp(1i*omega*(x(:,1) + x(:,2)));
    
    a = 1/2;
    b = 2*a;
    [node,elem] = squaremesh([-a,a,-a,a],h);
    [lnode,lelem] = squaremesh([-b,b,-b,b],h);
    if omega > 100*pi
        [node,elem] = squaremesh([-a,a,-a,a],4*h);
    end
    if omega > 400*pi
        [node,elem] = squaremesh([-a,a,-a,a],8*h);
    end 
    if omega > 800*pi
        [node,elem] = squaremesh([-a,a,-a,a],16*h);
    end 
    N = size(node,1);
    err1 = zeros(N,Nray);
    err2 = err1;
    rec_z = err1;
    exray = [pi/4*ones(N,1), 3*pi/4*ones(N,1),...
        5*pi/4*ones(N,1), 7*pi/4*ones(N,1)];
    
    tic;
    for i = 1:1 %N
        x0 = node(i,1);
        y0 = node(i,2);
        r = -wl:h:wl;
        r = r';
        est_ang = pi/3;
        x = r*cos(est_ang) + x0;
        y = r*sin(est_ang) + y0;
        xy = [x, y];
        u = exu(xy);
        
        %% MPM
        [z] = Matrix_Pencil_2(u);
        [~,sz] = sort(imag(z));
        z = z(sz);  % [5/4*pi, 7/4*pi, 3/4*pi, 1/4*pi]
        err1(i,:) = [abs(5*pi/4 - est_ang - acos(log(z(1))/(1i*omega*h)) ),...
            abs(2*pi - (7*pi/4 - est_ang) - acos(log(z(2))/(1i*omega*h)) ),...
            abs(3*pi/4 - est_ang - acos(log(z(3))/(1i*omega*h)) ),...
            abs(est_ang - pi/4 - acos(log(z(4))/(1i*omega*h)) )];
        
        %% NMLA
        c0 = 1;
        Rest = 4;
        plt = 0;
        uu = exu(lnode);
        ux = Ux(lnode);   uy = Uy(lnode); 
        [angle] = NMLA_2D_2nd(x0,y0,c0,omega,Rest,lnode,lelem,uu,ux,uy,[],1/8,Nray,'num',0,plt);
        err2(i,:) = angle - [1/4*pi, 3/4*pi, 5/4*pi, 7/4*pi];
%         err2(i) = abs(exp(1i*omega*h*cos(pi/4 - est_ang)) - rec_z(i));
    end
    toc;
    
    rec_err1(j) = norm(err1(:))*h/(norm(exray(:))*h);    
    rec_err2(j) = norm(err2(:))*h/(norm(exray(:))*h);
end

rec_err1
      
rec_err2

% showrate(rec_omega, rec_err1');


% acos(asin(max(imag(z)))/(omega*T))