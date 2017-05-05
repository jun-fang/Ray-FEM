function uint = interpolation2(x, y, u, xy)
%% interpolate u at xy nodes: valid for retangular mesh
% complexity O(ru x cu), where ru and cu are row and column numbers of u

xmin = min(x);    xmax = max(x);
ymin = min(y);    ymax = max(y);
% if min(xy(:,1))<xmin-eps || max(xy(:,1))>xmax+eps...
%         || min(xy(:,2))<ymin-eps || max(xy(:,2))>ymax+eps
%     fprintf('Attention!!! xy nodes out of the interpolation domain !!!\n');
%     return;
% end


cu = size(u,2);              % column number of u
rxy = size(xy,1);            % row number of xy
uint = zeros(rxy,cu);        % output

m = length(x);  n = length(y);
h = (xmax - xmin) / (m - 1);

ix = 1 + floor((xy(:,1)-xmin)/h + eps);     % find the location in x coordinate
iy = 1 + floor((xy(:,2)-ymin)/h + eps);     % find the location in y coordinate

ix = min(ix, m-1);  ix = max(ix, 1);
iy = min(iy, n-1);  iy = max(iy, 1);

% 4--3
%  | \  |
% 1--2
x = x(:);  y = y(:);
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
for i = 1: cu
    uint1 = a12.*findu(ix, iy+1, n, u(:,i)) + a24.*findu(ix, iy, n, u(:,i)) + a41.*findu(ix+1, iy, n, u(:,i)) ;
    uint2 = a23.*findu(ix, iy+1, n, u(:,i)) + a24.*findu(ix+1, iy+1, n, u(:,i)) + a34.*findu(ix+1, iy, n, u(:,i)) ;
    uint(:,i) = uint1.*( sum(xy,2) <= sum(xy1,2)+h) + uint2.*( sum(xy,2) > sum(xy1,2)+h) ;
end

end


function uxy = findu(ix, iy, n, u)
ind = (ix - 1)*n + iy;
uxy = u(ind,:);
end




