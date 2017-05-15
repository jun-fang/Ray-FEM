function [] = ray_field(ray,node,gap_num,pct)
%% plot the ray field 
% gap_num: gap number to plot the ray directions 
xmin = node(1,1);    xmax = node(end,1);
ymin = node(1,2);    ymax = node(end,2);
h = node(2,2) - node(1,2);    % mesh size
h = 1/round(1/h);
m = round((xmax-xmin)/h + 1);
n = round((ymax-ymin)/h + 1);
% pct = 1/10; rescale the length of the ray direction
if ((m-1)/gap_num)>floor(((m-1)/gap_num))...
        || ((n-1)/gap_num)>floor(((n-1)/gap_num))
    fprintf('Gap number is not correct!!!\n');
    return;
end

[X,Y] = meshgrid(xmin:h:xmax,ymin:h:ymax);
ch = gap_num*h;

dX = (X - xmin)/ch - round((X - xmin)/ch);
dY = (Y - ymin)/ch - round((Y - ymin)/ch);
dX = abs(dX(:));  dY = abs(dY(:));
cn = find((dX<h/10).*(dY<h/10));

if iscell(ray)    %% ray is a cell
    x = [];
    y = [];
    rayx = [];
    rayy = [];
        
    for i = 1:length(cn)
        ii = cn(i);
        ray_num = size(ray{ii},2);
        v = ones(ray_num,1);
        x = [x;node(ii,1)*v];
        y = [y;node(ii,2)*v];
        re = pct*real(ray{ii});
        im = pct*imag(ray{ii});
        rayx = [rayx; re(:)];
        rayy = [rayy; im(:)];
    end
%     quiver(x,y,rayx,rayy);    
else              %% ray is an array
    Nray = size(ray,2);
    x = node(cn,1);   y = node(cn,2);
    x = repmat(x,Nray,1);
    y = repmat(y,Nray,1);
    x = x(:);  y = y(:);
    rayx = pct*real(ray(cn,:));
    rayx = rayx(:);
    rayy = pct*imag(ray(cn,:));
    rayy = rayy(:);
%     quiver(x,y,rayx,rayy);
end
% d = abs(rayx) + abs(rayy);
% rayx = rayx + 1/100*(d<10*eps);
% rayy = rayy + 1/100*(d<10*eps);
quiver(x,y,rayx,rayy);
d = abs(rayx) + abs(rayy);
dx = x(d<10*eps);  dy = y(d<10*eps);
hold on; scatter(dx,dy,'b.')
title('Ray field','FontSize', 16);