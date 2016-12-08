
%% basic_cutoff
% x = -100:200;
% x = x/100;
% y = basic_cutoff(x);
% plot(x,y)


%% basic_cutoff_first_derivative
% x = -100:200;
% x = x/100;
% y = basic_cutoff_first_derivative(x);
% plot(x,y)

% a = 2*rand(1); b = 2*rand(1);
% int1 = basic_cutoff(b) - basic_cutoff(a);
% int2 = integral(@basic_cutoff_first_derivative,a,b);
% int1-int2
% 

%% basic_cutoff_second_derivative
% x = -100:200;
% x = x/100;
% y = basic_cutoff_second_derivative(x);
% plot(x,y)

% a = rand(1); b = rand(1);
% int1 = basic_cutoff_first_derivative(b) - basic_cutoff_first_derivative(a);
% int2 = integral(@basic_cutoff_second_derivative,a,b);
% a,b,int1-int2


%% Test the convergence order for cut-off functions wrt epsilon
% cut-off L2: \sim eps^1, L^{\infy} \sim eps^0
% gradient L2: \sim eps^0, L^{\infty} \sim eps^-1
% laplacian L2: \sim eps^-1, L^{\infty} \sim eps^-2

% h = 1/1900;  a = 1/2;
% [node,elem] = squaremesh([-a,a,-a,a],h);
% p = node;
% n = [1,2,3,4,5];
% x = 1.7.^-(n+2);
% y0 = x; y1 = x; y2 = x;
% for ni = 1:length(n)
%     ni
%     epsilon = x(ni);
%     a = epsilon;  b = 2*epsilon;
%     cg = cutoff_gradient(a,b,p,0,0);
%     abscg = abs(cg(:,1));%sqrt(abs(cg(:,1)).^2 + abs(cg(:,1)).^2);
%     cl = cutoff_laplacian(a,b,p,0,0);
%     cf = cutoff(a,b,p,0,0);
%     y0(ni) = norm(cf);
%     y1(ni) = norm(abscg);
%     y2(ni) = norm(cl);
% %     y0(ni) = norm(cf,inf);
% %     y1(ni) = norm(abscg,inf);
% %     y2(ni) = norm(cl,inf);
% end
% figure(2);
% subplot(3,1,1);
% show_convergence_rate(x,y0)
% subplot(3,1,2);
% show_convergence_rate(x,y1)
% subplot(3,1,3);
% show_convergence_rate(x,y2)



%% Test for convergence order for the right hand side
% fix epsilon, RHS \sim \omega^0.5, the first term dominates
%
% fix omega, the support of cut-off function depends on epsilon.   
%     RHS:  L^2 \sim eps^-0.5, L^{\infty} \sim eps^-1.5
%     RHS1: L^2 \sim eps^-0.5, L^{\infty} \sim eps^-1.5
%     RHS2: L^2 \sim eps^-1.5, L^{\infty} \sim eps^-2.5
% 
% when omega*eps^2 >= 18, the first term RHS1 dominates

h = 1/900;  a = 1/2;
xs = 0;  ys = 0;
[node,elem] = squaremesh([-a,a,-a,a],h); p = node;
xx = p(:,1)-xs;   yy = p(:,2)-ys;
r = sqrt(xx.^2 + yy.^2);  
% x = [80 160 320 640]*pi;

% x = [100 200 400 800]*pi;

omega = 640*pi;
n = [1,2,3,4,5];
x = 1.7.^-(n+2);

y0 = x;  y1 = x; y2 = x; 

% epsilon = 0.237;
% a = epsilon;  b = 2*epsilon;

for ni = 1:length(x)
    ni
    
    %     omega = x(ni);

    epsilon = x(ni); %sqrt(18/x(ni));
    a = epsilon;  b = 2*epsilon;


    ub = 1i/4*besselh(0,1,omega*r);
    ub_g1 = -1i/4*omega*besselh(1,1,omega*r)./r.*xx;  % partial derivative wrt x
    ub_g2 = -1i/4*omega*besselh(1,1,omega*r)./r.*yy;  % partial derivative wrt y

    cf_lap = cutoff_laplacian(a,b,p,xs,ys);
    rhs2 = ub.*cf_lap;
    rhs2(r<=a) = 0; rhs2(r>=b) = 0;
    
    cf_grad = cutoff_gradient(a,b,p,xs,ys);
    rhs1 = ub_g1.*cf_grad(:,1) + ub_g2.*cf_grad(:,2);
    rhs1(r<=a) = 0; rhs1(r>=b) = 0;
    
    rhs = rhs1 + rhs2;
    
    y0(ni) = norm(rhs,inf);
    y1(ni) = norm(rhs1,inf);
    y2(ni) = norm(rhs2,inf);
     
%     y0(ni) = norm(rhs)*h;
%     y1(ni) = norm(rhs1)*h;
%     y2(ni) = norm(rhs2)*h;
end
figure(3);
subplot(3,1,1);
show_convergence_rate(x,y0)
subplot(3,1,2);
show_convergence_rate(x,y1)
subplot(3,1,3);
show_convergence_rate(x,y2)




%% Test for the exact solution of smooth part
% h = 1/900;  a = 1/2;
% xs = 0;  ys = 0;
% [node,elem] = squaremesh([-a,a,-a,a],h); p = node;
% xx = p(:,1)-xs;   yy = p(:,2)-ys;
% r = sqrt(xx.^2 + yy.^2);
% 
% nt = 4;
% y = zeros(1,nt);
% x = [80 160 320 640]*pi;
% 
% for ii = 1:nt
%     
%     omega = x(ii);
%     ub = 1i/4*besselh(0,1,omega*r);
%     epsilon = 0.143;
%     a = epsilon;  b = 2*epsilon;
%     cf = cutoff(a,b,p,xs,ys);
%     
%     uex = (1-cf).*ub;
%     uex(r<=a) = 0;
%     
%     amp = uex./exp(1i*omega*r);
%     amp(r<=a) = 0;
%     l2amp = norm(amp)*h;
%     
%     y(ii) = l2amp;
% end
% 
% show_convergence_rate(x,y)


