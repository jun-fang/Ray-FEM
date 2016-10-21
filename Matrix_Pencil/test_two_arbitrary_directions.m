%% test for two arbitrary sampling directions
% four point source problem

global omega xs ys a;
omega = 50*pi;
xs = 100;  ys = 100;
a = 1/2;

% omega = sqrt(omega);
pde = MPM_data5;
Nray = 4;
NPW = 6;
h = 1/round((NPW*omega)/(2*pi));   % mesh size
wl = 2*pi/omega;  

[node,elem] = squaremesh([-a,a,-a,a],h);
N = size(node,1);
rec_err = zeros(N,Nray);

tic;
for i = 1:N

    i
    N
% d = 1;%200*(rand(1)-0.5);
% x0 = d*(rand(1)-0.5);  y0 = d*(rand(1)-0.5);
x0 = node(i,1); y0 =node(i,2);
exray = pde.ray_ang([x0,y0]);
r = -wl:h:wl;  r = r';

est_ang =  pi/6;% + (-1)^(exray(i)>pi)*pi/30;
x = r*cos(est_ang) + x0;  y = r*sin(est_ang) + y0;
xy = [x, y];
u = pde.ex_u(xy);
u = u.*sqrt((x-xs).^2 + (y-ys).^2).^(1/2);% + 10^(-3)*ones(size(u));
[z] = Matrix_Pencil(u,4);

dang1 = acos(1i*imag(log(z))/(1i*omega*h));
ang1 = [est_ang + dang1; est_ang - dang1];

ang1 = principal_angle(ang1);
% err1 = min(abs(abs(exray - est_ang) - abs(acos(1i*imag(log(z))/(1i*omega*h))) ));
err1 = min(abs(cos(exray' - est_ang) - 1i*imag(log(z))/(1i*omega*h)));

% fprintf('\n-----------------------------------\n')

est_ang =  pi/12;% + (-1)^(exray(i)>pi)*pi/30;
x = r*cos(est_ang) + x0;  y = r*sin(est_ang) + y0;
xy = [x, y];
u = pde.ex_u(xy);
u = u.*sqrt((x-xs).^2 + (y-ys).^2).^(1/2);
tol = 3;
L = 10;
[z] = Matrix_Pencil(u, 4);

dang2 = acos(1i*imag(log(z))/(1i*omega*h));
ang2 = [est_ang + dang2; est_ang - dang2];
ang2 = principal_angle(ang2);
% err1 = min(abs(abs(exray - est_ang) - abs(acos(1i*imag(log(z))/(1i*omega*h))) ));
err2 = min(abs(cos(exray' - est_ang) - 1i*imag(log(z))/(1i*omega*h)));
% d
% x0
% y0
% exray/pi
% min(abs(ang1 - exray))
% min(abs(ang2 - exray))
angle = find_angle(ang1, ang2);
rec_err(i,:) = abs(angle - exray);
end
toc;
norm(rec_err,inf)
        