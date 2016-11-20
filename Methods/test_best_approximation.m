global omega;
pde = Helmholtz_data1;
NPW = 6;
a = 1/2;

omegas = [10 20 40 80]*pi;
errors = 0*omegas;
tic
for i = 1:length(omegas)
    i
    tic;
    omega = omegas(i);
    h = 1/(NPW*round(omega/(2*pi)));
    [node,elem] = squaremesh([-a,a,-a,a],h);
    ray = ex_ray(node,2,2,1);
    [error] = best_approximation(node,elem,ray,omega,pde);
    errors(i) = error;
    toc;
end
toc;

showrate(omegas,errors)