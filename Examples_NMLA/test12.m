function y = test12(x,s,x0,y0)
y = y0.*(x./x0).^s;
y = mean(y);

if (0)
global omega xs ys;

%% Test 1
pde1 = Helmholtz_data1;
xs = 2; ys = 2;
cp_omega1 = [40,48,64,80,104,128,160,200,240,280,320,360]*pi;
% cp_omega1 = cp_omega1;
NPW = 6;
exu_norm1 = 0*cp_omega1;
ray_norm1 = 0*cp_omega1;
for i = 1:length(cp_omega1)
    omega = cp_omega1(i);
    h = 1/(NPW*round(omega/(2*pi)));
    [node,elem] = squaremesh([-1/2,1/2,-1/2,1/2],h);
    exu = pde1.ex_u(node);
    exu_norm1(i) = h*norm(exu,2);
    exray = ex_ray_angle(node,xs,ys);
    ray_norm1(i) = (h*norm(exray,2));
end
idx1 = [1 4 7 11];
idx2 = [1 6 16 36];
load('result1_wrt_omega.mat');
omega = rec_omega1(idx2);
h = rec_h(idx2);
rayerr1 = rec_ang_err1(idx2).*ray_norm1(idx1);
rayerr2 = rec_ang_err2(idx2).*ray_norm1(idx1);
uerr_ray1 = rec_NR_err1(idx2).*exu_norm1(idx1);
uerr_ray2 = rec_NR_err2(idx2).*exu_norm1(idx1);
uerr_exray = rec_ER_err(idx2).*exu_norm1(idx1);
save('result1_errs.mat','omega','h','rayerr1','rayerr2','uerr_ray1','uerr_ray2','uerr_exray');
fprintf( '\n\nAngle L2 error 1:       ');
fprintf( '&  %1.2d  ',rayerr1);
fprintf( '\n\nAngle L2 error 2:       ');
fprintf( '&  %1.2d  ',rayerr2);
fprintf( '\n\nNR-FEM L2 error 1:       ');
fprintf( '&  %1.2d  ',uerr_ray1);
fprintf( '\n\nNR-FEM L2 error 2:       ');
fprintf( '&  %1.2d  ',uerr_ray2);
fprintf( '\n\nER-FEM L2 error 1:       ');
fprintf( '&  %1.2d  ',uerr_exray);


%% Test 2
pde2 = Helmholtz_data2;
xs = 20; ys = 20;
cp_omega2 = [10 20 30 40 60 80 100]*pi;
cp_omega2 = cp_omega2';
NPW = 6;
exu_norm2 = 0*cp_omega2;
ray_norm2 = 0*cp_omega2;
for i = 1:7
    omega = cp_omega2(i);
    h = 1/(NPW*round(omega/(2*pi)));
    [node,elem] = squaremesh([-1/2,1/2,-1/2,1/2],h);
    exu = pde2.ex_u(node);
    exu_norm2(i) = h*norm(exu,2);
    ray = pde2.ray(node);
    ray = ray(:);
    ray = [real(ray), imag(ray)];
    ray_dir = atan2(ray(:,2),ray(:,1));
    neg_index = find(ray_dir<0);
    ray_dir(neg_index) = ray_dir(neg_index) + 2*pi;
    ray_norm2(i) = (h*norm(ray_dir,2));
end
end