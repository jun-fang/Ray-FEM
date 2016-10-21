%% test 4: MPM for one point source problem, outside domain
% with numerically computed data
% low frequency


% rec_omega = [ 3200 6400]*pi;
rec_omega = [25 50 100 200 400 800 1600]*pi;
rec_err1 = zeros(length(rec_omega),1);
rec_err2 = rec_err1;
rec_err3 = rec_err1;
rec_err0 = rec_err1;
rec_err4 = rec_err1;

global omega;
global xs;
global ys;
global a;

pde = Helmholtz_data1_0;
plt = 0;                           % plot solution or not
fquadorder = 3;                    % numerical quadrature order
solver = 'DIR';

xs = 10;
ys = 10;

nmla = 0;
nn = 5;

for j = 1:nn%length(rec_omega)
    tic;
    j
    low_omega = sqrt(rec_omega(j));
    omega = rec_omega(j);
    NPW = 8;
    h = 1/round((NPW*omega)/(2*pi));   % mesh size
    
    
%     pde = MPM_data4;
    Nray = 1;
    
    a = 1/2;
    
    for ii = 1: 9
        if omega > 50*pi*2^(ii-1)
            [node,elem] = squaremesh([-a,a,-a,a],2^(ii)*2*h);
            nh = 2^(ii)*2*h;
        end
    end
    
    if omega <= 50*pi
        [node,elem] = squaremesh([-a,a,-a,a],4*h);
        nh = 4*h;
    end
    
    omega = low_omega;
    wl = 2*pi/omega;    % wavelength
    rh = wl/NPW;
    
    a = a + h*(round(wl/h) + 1);
    [lnode,lelem] = squaremesh([-a,a,-a,a],h);
    [lu,~,~,rel_L2_err] = Standard_FEM_IBC(lnode,lelem,omega,pde,fquadorder,solver,plt);
%     rel_L2_err
     eu = pde.ex_u(lnode);

    rec_err0(j) = norm(lu-eu,inf);

%         lu = pde.ex_u(lnode);
    
    %     [X,Y] = meshgrid(-b:h:b,-b:h:b);
    %     V = pde.ex_u([X(:),Y(:)]);
    %     V = reshape(V,size(X));
    %     V = V/sqrt(omega);
    
    
    
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
        d0 = sqrt((x0-xs)^2 + (y0-ys)^2);
        if d0 > 0*wl
            r = -wl:rh:wl;
            r = r';
            %             est_ang =  exray(i)/1.1 + pi/12;
            est_ang =  exray(i);% + (-1)^(exray(i)>pi)*pi/30;
            x = r*cos(est_ang) + x0;
            y = r*sin(est_ang) + y0;
            xy = [x, y];
            %             u = pde.ex_u(xy);
            %             nu = interpolation(lnode,lelem,xy,lu);
            eu = pde.ex_u(xy);
            
            nu = interpolation(lnode,lelem,xy,lu);
            %             tic;
            %             nu = interp2(X,Y,V,x,y,'cubic');
            %             toc;
            
            err3(i,:) = norm(nu-eu,inf);
            u = nu;
            u = u.*sqrt((x-xs).^2 + (y-ys).^2).^(1/2);
            
            %% MPM
            [z] = Matrix_Pencil_2(u);
            %             [~,sz] = sort(imag(z));
            %             z = z(sz);  % [5/4*pi, 7/4*pi, 3/4*pi, 1/4*pi]
            %             z = z(end);
            
            err1(i,:) = min(abs(abs(exray(i) - est_ang) - abs(acos(1i*imag(log(z))/(1i*omega*rh))) ));
            err2(i,:) = min(abs(cos(exray(i) - est_ang) - 1i*imag(log(z))/(1i*omega*rh)));
            
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
    
    
    
    rec_err1(j) = norm(err1(:))/norm(exray(:));
    rec_err2(j) = norm(err1(:),inf);
    rec_err3(j) = norm(err2(:))*nh;
    rec_err4(j) = norm(err3(:),inf);
    
    toc;
end

fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('xs = %d,    ys = %d   \n\n', xs,ys);
fprintf('MPM error at ');
fprintf('Sample data angle: %d pi ', est_ang/pi);
fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('Error_0     ');
fprintf('  &  %.2e',rec_err0');
fprintf('\n\n');

fprintf('Error_1     ');
fprintf('  &  %.2e',rec_err1');
fprintf('\n\n');

fprintf('Error_2     ');
fprintf('  &  %.2e',rec_err2');
fprintf('\n\n');

fprintf('Error_3     ');
fprintf('  &  %.2e',rec_err3');
fprintf('\n\n');

fprintf('Error_4     ');
fprintf('  &  %.2e',rec_err4');
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

subplot(1,2,1);
showrate(rec_omega(1:nn),rec_err1(1:nn)');
xlabel('\omega');
ylabel('Err_{\theta, {L^{2}(\Omega)}}');


subplot(1,2,2);
showrate(rec_omega(1:nn),rec_err2(1:nn)');
xlabel('\omega');
ylabel('Err_{\theta, {L^{\infty}(\Omega)}}');





% figure(1);
% showrate(rec_omega, rec_err1');
%
% figure(2);
% showrate(rec_omega, rec_err2');
%
% figure(3);
% showrate(rec_omega, rec_err3');