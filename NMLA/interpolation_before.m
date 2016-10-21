function uint = interpolation_before(node,elem,xy,u)
%% interpolate u at xy nodes

cu = size(u,2);              % column number of u
rxy = size(xy,1);            % row number of xy
uint = zeros(rxy,cu);        % output
N = size(node,1);            % number of nodes
NT = size(elem,1);           % number of elements
n = round(sqrt(N));          % number of grid points in each dimension
h = node(2,2) - node(1,2);   % grid width 

xmin = node(1,1);   ymin = node(1,2);
xyn = xy;                    % record the index of x,y coordinates
xyn(:,1) = floor((xy(:,1)-xmin)/h + eps);
xyn(:,2) = floor((xy(:,2)-ymin)/h + eps);

xyn = min(xyn,n-2);
ind = xyn(:,1)*(n - 1) + xyn(:,2) + 1;
temp = xy - xyn*h;
ind = ind + (temp(:,1)<temp(:,2)-eps)*NT/2;

ceara = h^2/2;

max(ind)
a1 = tri_area(node(elem(ind,2),:),node(elem(ind,3),:),xy)/ceara;
a2 = tri_area(node(elem(ind,3),:),node(elem(ind,1),:),xy)/ceara;
a3 = tri_area(node(elem(ind,1),:),node(elem(ind,2),:),xy)/ceara;



for i = 1:cu
    uint(:,i) = a1.*u(elem(ind,1),i) + a2.*u(elem(ind,2),i) + a3.*u(elem(ind,3),i);
end

end