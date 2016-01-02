function uray = ray_solution(node,elem,omega,speed,v,ray,xy)
Nray = size(ray,2);

xmin = node(1,1);    xmax = node(end,1);
ymin = node(1,2);    ymax = node(end,2);
if min(xy(:,1))<xmin-eps || max(xy(:,1))>xmax+eps...
        || min(xy(:,2))<ymin-eps || max(xy(:,2))>ymax+eps
    fprintf('Attention!!! xy nodes out of the interpolation domain !!!\n');
    return;
end

rxy = size(xy,1);            % row number of xy
uray = zeros(rxy,1);

N = size(node,1);            % number of nodes
NT = size(elem,1);           % number of elements
n = round(sqrt(N));          % number of grid points in each dimension
h = node(2,2) - node(1,2);   % grid width 
each_area = h^2/2;           % area of each triangle
v = reshape(v,N,Nray);


xyn = xy;                    % record the index of x,y coordinates
x = xy(:,1);
y = xy(:,2);
xyn(:,1) = floor((x-xmin)/h + eps);     % find the location in x coordinate
xyn(:,2) = floor((y-ymin)/h + eps);     % find the location in y coordinate

xyn = min(xyn,n-2);        
ind = xyn(:,1)*(n - 1) + xyn(:,2) + 1;        % find the element index where xy locates
temp = xy - xyn*h;
ind = ind + (temp(:,1)<temp(:,2)-eps)*NT/2;   % modify the index

% compute the component of each little triangle
n1 = elem(ind,1); n2 = elem(ind,2); n3 = elem(ind,3); 
node1 = node(n1,:); node2 = node(n2,:); node3 = node(n3,:);
a1 = tri_area(node2,node3,xy)/each_area;
a2 = tri_area(node3,node1,xy)/each_area;
a3 = tri_area(node1,node2,xy)/each_area;

for i = 1:Nray
    re1 = real(ray(n1,i)); im1 = imag(ray(n1,i));
    re2 = real(ray(n2,i)); im2 = imag(ray(n2,i));
    re3 = real(ray(n3,i)); im3 = imag(ray(n3,i));
    exp1 = exp(1i*omega./speed(node1).* ( re1.*x + im1.*y ));
    exp2 = exp(1i*omega./speed(node2).* ( re2.*x + im2.*y ));
    exp3 = exp(1i*omega./speed(node3).* ( re3.*x + im3.*y ));
    uray = uray + v(n1,i).*a1.*exp1 + v(n2,i).*a2.*exp2 + v(n3,i).*a3.*exp3;
end



