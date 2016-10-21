function uint = interpolation(node,elem,xy,u)
%% interpolate u at xy nodes
% complexity O(ru x cu), where ru and cu are row and column numbers of u

xmin = node(1,1);    xmax = node(end,1);
ymin = node(1,2);    ymax = node(end,2);
if min(xy(:,1))<xmin-eps || max(xy(:,1))>xmax+eps...
        || min(xy(:,2))<ymin-eps || max(xy(:,2))>ymax+eps
    fprintf('Attention!!! xy nodes out of the interpolation domain !!!\n');
    return;
end

cu = size(u,2);              % column number of u
rxy = size(xy,1);            % row number of xy
uint = zeros(rxy,cu);        % output
N = size(node,1);            % number of nodes
NT = size(elem,1);           % number of elements
n = round(sqrt(N));          % number of grid points in each dimension
h = node(2,2) - node(1,2);   % grid width 
each_area = h^2/2;           % area of each triangle

xyn = xy;                    % record the index of x,y coordinates
xyn(:,1) = floor((xy(:,1)-xmin)/h + eps);     % find the location in x coordinate
xyn(:,2) = floor((xy(:,2)-ymin)/h + eps);     % find the location in y coordinate

xyn = min(xyn,n-2);        
ind = xyn(:,1)*(n - 1) + xyn(:,2) + 1;        % find the element index where xy locates
temp = xy - xyn*h;
ind = ind + (temp(:,1)<temp(:,2)-eps)*NT/2;   % modify the index

% compute the component of each little triangle
a1 = tri_area(node(elem(ind,2),:),node(elem(ind,3),:),xy)/each_area;
a2 = tri_area(node(elem(ind,3),:),node(elem(ind,1),:),xy)/each_area;
a3 = tri_area(node(elem(ind,1),:),node(elem(ind,2),:),xy)/each_area;

% interpolation
for i = 1:cu
    uint(:,i) = a1.*u(elem(ind,1),i) + a2.*u(elem(ind,2),i) + a3.*u(elem(ind,3),i);
end

end