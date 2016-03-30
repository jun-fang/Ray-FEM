function ray = ex_ray(node,xs,ys,opt)
%% exact ray angle for point source in homogeneous medium
xx = node(:,1)-xs; 
yy = node(:,2)-ys;
ray = atan2(yy,xx);
d = sqrt(xx.^2 + yy.^2);
ray = ray + 2*pi*(ray < 0);
if opt 
    ray = exp(1i*ray).*(d>10*eps);
end
