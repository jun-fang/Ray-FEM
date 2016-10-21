function f = approximate_delta_2(xc,yc,h,xy)
x = xy(:,1); y = xy(:,2);
n = length(x);
idx1 = (x-xc).*(y-yc)>0;
mx = max(abs(x-xc), abs(y-yc));
idx2 = mx<h;
node_north = [xc*ones(n,1),(yc+h)*ones(n,1)];  
node_south = [xc*ones(n,1),(yc-h)*ones(n,1)];
node_west = [(xc-h)*ones(n,1),yc*ones(n,1)];
node_east = [(xc+h)*ones(n,1),yc*ones(n,1)];
area1 = tri_area(xy,node_west,node_north);
area2 = tri_area(xy,node_east,node_south);
m1 = min(area1,area2)/(1/2*h*h);
m2 = max(area1,area2)/(1/2*h*h);
f = (1 - mx/h).*idx1.*idx2 + m1.*(m1<=1).*(m2<=2).*(1-idx1);