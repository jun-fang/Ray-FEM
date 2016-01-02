function ray = ex_ray(node,xs,ys)
%% exact ray angle for point source in homogeneous medium
xx = node(:,1)-xs; 
yy = node(:,2)-ys;
ray = atan2(yy,xx);
ray = ray + 2*pi*(ray < 0);