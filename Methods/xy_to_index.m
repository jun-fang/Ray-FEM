%% map from xy coordinates to the index of grid nodes

function [index, x, y] = xy_to_index(node,x0,y0)

xmin = node(1,1);    % xmax = node(end,1);
ymin = node(1,2);    ymax = node(end,2);

h = node(2,2) - node(1,2);    % mesh size
n = round((ymax-ymin)/h + 1);

xn = round((x0 - xmin)/h);   
yn = round((y0 - ymin)/h);
index = n*xn + yn + 1;          
x = xmin + xn*h;
y = ymin + yn*h;
% h = node(2,2) - node(1,2);
% d = sqrt((node(:,1)-x0).^2 + (node(:,2)-y0).^2 );
% index = find(d < h/2);
end

