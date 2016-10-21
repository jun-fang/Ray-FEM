%% Test for numerical quadrature in triangle with vertices 
% (x0,y0) (x0+h,y0) (x0,y0+h)
    
omega = 40*pi;
NPW = 6;                           % number of grid points per wavelength
h = 1/round((NPW*omega)/(2*pi));
x0 = rand(1); 
y0 = rand(1);
xy0 = [x0,y0]
ang0 = 5*pi/4;
d0 = [cos(ang0),sin(ang0)];

ang0 = ang0 -pi/2;
d1 = [cos(ang0),sin(ang0)];

d0 = d0 - d1;

f = @(x) (  exp(1i*omega*( dot(d0,x-xy0) ))    );

xy1 = [x0,y0];
xy2 = [x0+h,y0];
xy3 = [x0,y0+h];
area = h^2/2;

%% 3rd Numerical Quadrature
fquadorder = 2;  
[lambda,weight] = quadpts(fquadorder);
phi = lambda;           % the value of piesewise linear basis functions at numerical quadrature nodes
nQuad = size(lambda,1);
I = 0;
for p = 1:nQuad
    pxy = lambda(p,1)*xy1 ...
        + lambda(p,2)*xy2 ...
        + lambda(p,3)*xy3;
    (pxy - xy0)/h
    fp = f(pxy);
    I = I + weight(p)*fp;
end
itg1 = I*area;

fquadorder = 9;  
[lambda,weight] = quadpts(fquadorder);
phi = lambda;           % the value of piesewise linear basis functions at numerical quadrature nodes
nQuad = size(lambda,1);
I = 0;
for p = 1:nQuad
    pxy = lambda(p,1)*xy1 ...
        + lambda(p,2)*xy2 ...
        + lambda(p,3)*xy3;
    fp = f(pxy);
    I = I + weight(p)*fp;
end
itg2 = I*area;


omega^2*(itg1 - itg2)
















