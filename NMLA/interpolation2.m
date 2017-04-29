function uint = interpolation2(x, y, u, xy)
%% interpolate u at xy nodes
% complexity O(ru x cu), where ru and cu are row and column numbers of u

xmin = min(x);    xmax = max(x);
ymin = min(y);    ymax = max(y);
if min(xy(:,1))<xmin-eps || max(xy(:,1))>xmax+eps...
        || min(xy(:,2))<ymin-eps || max(xy(:,2))>ymax+eps
    fprintf('Attention!!! xy nodes out of the interpolation domain !!!\n');
    return;
end

n = length(x);
h = (xmax - xmin) / (n - 1);

ix = 1 + floor((xy(:,1)-xmin)/h + eps);     % find the location in x coordinate
iy = 1 + floor((xy(:,2)-ymin)/h + eps);     % find the location in y coordinate

% 4--3
%  | \  |
% 1--2   

xy1 = [x(ix), y(iy)];
xy2 = [x(ix+1), y(iy)];
xy3 = [x(ix+1), y(iy+1)];
xy4 = [x(ix), y(iy+1)];

% compute the component of each little triangle
area = 1/2*h*h;
a12 = tri_area(xy1, xy2, xy)/area;
a23 = tri_area(xy2, xy3, xy)/area;
a24 = tri_area(xy2, xy4, xy)/area;
a34 = tri_area(xy3, xy4, xy)/area;
a41 = tri_area(xy4, xy1, xy)/area;

% interpolation
uint1 = a12.*u(iy+1, ix) + a24.*u(iy, ix) + a41.*u(iy, ix+1);
uint2 = a23.*u(iy+1, ix) + a24.*u(iy+1, ix+1) + a34.*u(iy, ix+1);
uint = uint1.*( sum(xy,2) <= sum(xy1)+h) + uint2.*( sum(xy,2) > sum(xy1)+h) ;


end