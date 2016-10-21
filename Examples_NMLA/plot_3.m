figure(26);

subplot(3,2,1);
mesh(theta/pi,rn,real(u6));
xlabel('\theta/\pi');
ylabel('r/\lambda');
% axis equal;
axis tight;
az = 0;
el = 90;
view(az, el);
title('Ray-FEM NPW = 6');
        
subplot(3,2,2);
mesh(theta/pi,rn,real(us6));
xlabel('\theta/\pi');
ylabel('r/\lambda');
% axis equal;
axis tight;
az = 0;
el = 90;
view(az, el);
title('S-FEM NPW = 6');

subplot(3,2,3);
mesh(theta/pi,rn,real(u8));
xlabel('\theta/\pi');
ylabel('r/\lambda');
% axis equal;
axis tight;
az = 0;
el = 90;
view(az, el);
title('Ray-FEM NPW = 8');

subplot(3,2,4);
mesh(theta/pi,rn,real(us8));
xlabel('\theta/\pi');
ylabel('r/\lambda');
% axis equal;
axis tight;
az = 0;
el = 90;
view(az, el);
title('S-FEM NPW = 8');

subplot(3,2,5);
mesh(theta/pi,rn,real(u10));
xlabel('\theta/\pi');
ylabel('r/\lambda');
% axis equal;
axis tight;
az = 0;
el = 90;
view(az, el);
title('Ray-FEM NPW = 10');

subplot(3,2,6);
mesh(theta/pi,rn,real(us10));
xlabel('\theta/\pi');
ylabel('r/\lambda');
% axis equal;
axis tight;
az = 0;
el = 90;
view(az, el);
title('S-FEM NPW = 10');