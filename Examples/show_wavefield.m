function [] = show_wavefield(a,h,u,y)
n = round(2*a/h) + 1;
u = reshape(u,n,n);
y = h*ceil((y+a)/h);
yn = round(y/h) + 1;
xx = -a:h:a;
uh = u(yn,:);

plot(xx,real(uh),'r');

xlabel('x');
ylabel('Real part wavefield');
% legend('Ray-FEM solution','S-FEM solution','Reference solution','LOCATION','Best');
title(['Wavefield at y = ' num2str(y-a)],'FontSize', 14)

