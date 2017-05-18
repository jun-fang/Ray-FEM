% domain
sdx = 1.5; sdy = 0.5;
small_domain = [-sdx, sdx, -sdy, sdy];
[node,elem] = squaremesh(small_domain,h);


xmin = min(node(:,1)); xmax = max(node(:,1));
ymin = min(node(:,2)); ymax = max(node(:,2));
h = node(2,2) - node(1,2);
m = round((xmax - xmin)/h) + 1;
n = round((ymax - ymin)/h) + 1;

rnode = 0*node;   
rnode(:,1) = (node(:,1) - xmin)*4000/1000;
rnode(:,2) = (ymax - node(:,2))*4000/1000;

rxmin = min(rnode(:,1)); rxmax = max(rnode(:,1));
rymin = min(rnode(:,2)); rymax = max(rnode(:,2));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot the wave speed 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% rspeed = speed(node)*1500;
% smin = min(rspeed); smax = max(rspeed);
% 
% figure('position', [100, 100, 1650, 600]);
% showsolution(rnode,elem,rspeed,2); colorbar; axis equal; axis tight;
% xlabel('x[km]', 'FontSize', 20);  ylabel('y[km]','FontSize',20);
% xlim([0.999999*rxmin 1.000001*rxmax]); 
% ylim([0.999999*rymin 1.000001*rymax]);
% caxis([0.999999*smin 1.000001*smax]);
% set(gca,'position',[0.05 0.125 0.875 0.85],'units','normalized')
% set(gca,'YDir','reverse')
% set(gca,'fontsize',18)
% 
% filename = 'ex4_Marmousi_wavespeed';
% print(filename,'-depsc','-r500');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot the wave field 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ru = real(uh2);
umin = min(ru); umax = max(ru);

figure('position', [100, 100, 1650, 600]);
showsolution(rnode,elem,ru,2); colorbar; axis equal; axis tight;
xlabel('x[km]', 'FontSize', 20);  ylabel('y[km]','FontSize',20);
xlim([0.999999*rxmin 1.000001*rxmax]); 
ylim([0.999999*rymin 1.000001*rymax]);
caxis([0.999999*umin 1.000001*umax]);
set(gca,'position',[0.05 0.125 0.875 0.85],'units','normalized')
set(gca,'YDir','reverse')
set(gca,'fontsize',18)

filename = 'ex4_Marmousi_wavefield_100pi_NPW_4';
print(filename,'-depsc','-r500');



