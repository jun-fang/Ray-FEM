function idx = Marmousi_index(node, h)
x0 = -1.75;  y0 = 0.75;
yn = round(1.5/h) + 1;
ix = round( (node(:,1) - x0) / h );
iy = round( (y0 - node(:,2)) / h );
idx = ix*yn + iy + 1;
