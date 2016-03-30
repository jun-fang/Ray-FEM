global omega;
global xs;
global ys;

xs = 0;
ys = 0;

omega = 3200*pi;
NPW = 6;
h = 1/round((NPW*omega)/(2*pi));   % mesh size
wl = 2*pi/omega;    % wavelength

pde = MPM_data4;

x0 = rand(1);   
y0 = rand(1);

r = -wl:h:wl;
r = r';
est_ang =  pde.ray([x0,y0]);
x = r*cos(est_ang) + x0;
y = r*sin(est_ang) + y0;
xy = [x, y];
u = pde.ex_u(xy);
u = u.*sqrt((x-xs).^2 + (y-ys).^2).^(1/2);

%% MPM
[z] = Matrix_Pencil_2(u);
[~,sz] = sort(imag(z));
z = z(sz);  % [5/4*pi, 7/4*pi, 3/4*pi, 1/4*pi]

abs(z-exp(1i*omega*h))