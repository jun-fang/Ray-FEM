a = 1/2;
h = 1/400;
[node,elem] = squaremesh([-a,a,-a,a],h);

xs = -0.4; ys = -0.4;
sigma = 0.15;


window = @(y,alpha, beta) 1*(abs(y) <= beta) + (abs(y) > beta).*(abs(y) < alpha)...
    .*exp(2*exp(-(alpha- beta)./(abs(y)-beta))./ ((abs(y)-beta)./(alpha- beta)-1) ) ;


xHet = 0.1;
yHet = 0.1;
nu = @(x,y) -0.5*exp( -1/(2*sigma^2)*((x-xHet).^2 + (y-yHet).^2) )...
    .*window(sqrt((x-xHet).^2 + (y-yHet).^2), 0.3,0.3  );

x = node(:,1);   y =  node(:,2);
nn = nu(x,y);

r = sqrt((x-xs).^2 + (y-ys).^2);
epsilon = 50/(80*pi);               % cut-off parameter

nn = nn + (r<= 2*epsilon);

showsolution(node,elem,nn,2);
colorbar;