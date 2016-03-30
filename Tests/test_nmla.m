pde = Helmholtz_data1;
x0 = 0; y0 = 0; c0 = 1;
global omega;
omega = 1000;
Rest = 4;
Nray = 1;
node = []; 
elem = [];
plt = 0;

h = 1/50;
a = 1;
b = 1;   
[node,elem] = squaremesh([-a,a,-b,b],h);
u = pde.ex_u(node);
speed = pde.speed;


[angs,Bs] = NMLA_2D(x0,y0,c0,omega,Rest,node,elem,0,0,0,pde,1/5,Nray,'ex',plt);

angs - pi/4
Bs - (exp(x0) + exp(y0))*exp(1i*omega*sqrt(2)/2*(x0+y0))  
Bs - pde.ex_u([x0,y0])

[opt_angs,opt_Bs,opt_res,niter] = gradient_descent(node,u,omega,speed,x0,y0,angs,Bs)
