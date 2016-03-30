function [] = ray_field(ray,node,gap_num)
%% plot the ray field 
% gap_num: gap number to plot the ray directions 
xmin = node(1,1);    xmax = node(end,1);
ymin = node(1,2);    ymax = node(end,2);
h = node(2,2) - node(1,2);    % mesh size
h = 1/round(1/h);
m = round((xmax-xmin)/h + 1);
n = round((ymax-ymin)/h + 1);

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
        re = 1/8*real(ray{ii});
        im = 1/8*imag(ray{ii});
        rayx = [rayx; re(:)];
        rayy = [rayy; im(:)];
    end
    quiver(x,y,rayx,rayy);    
else              %% ray is an array
    Nray = size(ray,2);
    x = node(cn,1);   y = node(cn,2);
    x = repmat(x,Nray,1);
    y = repmat(y,Nray,1);
    rayx = 1/8*real(ray(cn,:));
    rayx = rayx(:);
    rayy = 1/8*imag(ray(cn,:));
    rayy = rayy(:);
    quiver(x,y,rayx,rayy);
end
