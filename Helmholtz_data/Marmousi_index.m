function idx = Marmousi_index(node, x0, y0, h)
% make sure (x0, y0) is the upper right corner
x0 = -abs(x0);  y0 = abs(y0);
yn = round(2*y0/h) + 1;
ix = round( (node(:,1) - x0) / h );
iy = round( (y0 - node(:,2)) / h );
idx = ix*yn + iy + 1;
