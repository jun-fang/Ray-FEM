function [] = show_ray_solution(a,h,rh,u,su,ru,y)
n = round(2*a/h) + 1;
u = reshape(u,n,n);
su = reshape(su,n,n);


rn = round(2*a/rh) + 1;
ru = reshape(ru,rn,rn);

y = h*ceil((y+a)/h);
yn = round(y/h) + 1;
xx = -a:h:a;
uh = u(yn,:);
suh = su(yn,:);

ryn = round(y/rh) + 1;
rxx = -a:rh:a;
ruh = ru(ryn,:);

% -a + (yn-1)*h
% -a + (ryn-1)*rh


hold off;
plot(xx,real(uh),'r');
hold on;
plot(xx,real(suh),'b:');
hold on;
plot(rxx,real(ruh),'k--');
xlabel('x');
ylabel('Wavefield');
legend('Ray-FEM solution','S-FEM solution','Reference solution','LOCATION','Best');
title(['Wavefield at y = ' num2str(y-a)],'FontSize', 14)






% du = uh - uh2;
% rdu = ruh - ruh2;
% 
% norm(du,inf)
% norm(rdu,inf)






% uh2 = uh(end:-1:1);
% ruh2 = ruh(end:-1:1);






% figure(3);
% hold off;
% subplot(3,1,1)
% plot(xx,real(uh),'r-');
% hold on;
% plot(rxx,real(ruh),'k');
% xlabel('x');
% ylabel('Wavefield');
% legend('Ray-FEM solution','Reference solution','LOCATION','Best');
% title(['Wavefield at y = ' num2str(y-a)],'FontSize', 14)
% 
% subplot(3,1,2)
% plot(rxx,real(ruh),'k');
% hold on;
% 
% plot(rxx,real(ruh2),'r-');
% xlabel('x');
% ylabel('Wavefield');
% legend('Reference solution','Reflected reference solution','LOCATION','Best');
% title(['Wavefield at y = ' num2str(y-a)],'FontSize', 14)
% 
% 
% subplot(3,1,3)
% plot(xx,real(uh),'r-');
% hold on;
% 
% plot(xx,real(uh2),'k');
% xlabel('x');
% ylabel('Wavefield');
% legend('Ray-FEM solution','Reflected Ray-FEM solution','LOCATION','Best');
% title(['Wavefield at y = ' num2str(y-a)],'FontSize', 14)
