%% Plot
figure(20);   % wave-field
yy = 0.4;
yy = h*ceil((yy+a)/h);

% reference solution
rn = round(sqrt(length(ru)));
ref_u = reshape(ru,rn,rn);
ryn = round(yy/rh) + 1;
rxx = -a:rh:a;
ruh = ref_u(ryn,:);

plot(rxx,real(ruh),'k');
hold on;

% solution 1
n = round(sqrt(size(node,1)));
yn = round(yy/h) + 1;
xx = -a:h:a;
% uh = reshape(fu,n,n);
% uh = uh(yn,:);
% plot(xx,real(uh),'c^-');
% hold on;
uh = reshape(u,n,n);
uh = uh(yn,:);
plot(xx,real(uh),'ro-');
hold on;
uh = reshape(su,n,n);
uh = uh(yn,:);
plot(xx,real(uh),'b+:');
hold on;
   
xlabel('x');
ylabel('Real part wavefield');
legend('Reference solution','Ray-FEM solution with original RHS','Ray-FEM with more accurate RHS','S-FEM solution','LOCATION','Best');
title(['Wavefield at y = ' num2str(yy-a)],'FontSize', 14)
