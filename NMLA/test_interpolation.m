%% test for interpolation on retangular mesh

xs = 1; ys = 1;
a = 1/2; h = 1/500; h_c = 1/2;
domain1 = [-a,a,-a,a];
domain2 = [-2*a,2*a,-a,a];
[cnode, celem] = squaremesh(domain1,h_c);

figure(1);
showmesh(cnode,celem);
axis on
findnode(cnode);       % plot indices of all vertices
findelem(cnode,celem);  % plot indices of all triangles


[cnode, celem] = squaremesh(domain2,h_c);

figure(2);
showmesh(cnode,celem);
axis on
findnode(cnode);       % plot indices of all vertices
findelem(cnode,celem);  % plot indices of all triangles


[node, elem] = squaremesh(domain2,h);

h_c = 1/20;
[cnode, celem] = squaremesh(domain2,h_c);
cray = ex_ray(cnode,xs,ys,1);
ray = interpolation(cnode,celem,node,cray);
exray = ex_ray(node,xs,ys,1);
norm(ray-exray,inf)


h_c = 1/40;
[cnode, celem] = squaremesh(domain2,h_c);
cray = ex_ray(cnode,xs,ys,1);
ray = interpolation(cnode,celem,node,cray);
exray = ex_ray(node,xs,ys,1);
norm(ray-exray,inf)


h_c = 1/80;
[cnode, celem] = squaremesh(domain2,h_c);
cray = ex_ray(cnode,xs,ys,1);
ray = interpolation(cnode,celem,node,cray);
exray = ex_ray(node,xs,ys,1);
norm(ray-exray,inf)








