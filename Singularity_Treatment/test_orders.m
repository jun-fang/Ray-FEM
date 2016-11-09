%% Test for asymptotic order of hankel function 
% H(z) \sim O(z^-0.5)

% n = [1,2,3,4,5];
% x = 2.^n;
% y1 = besselh(0,1,x);
% y2 = besselh(1,1,x);
% showrate(x,abs(y1))
% showrate(x,abs(y2))


%% Test for cut-off function wrt epsilon
% gradient L2: \sim eps^-0.5, L^{ifty} \sim eps^-2
% laplacian L2: \sim eps^-2.5, L^{ifty} \sim eps^-4

% h = 1/3900;  a = 1/2;
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
%     abscg = sqrt(abs(cg(:,1)).^2 + abs(cg(:,1)).^2);
%     cl = cutoff_laplacian(a,b,p,0,0);
%     cf = cutoff(a,b,p,0,0);
%     y0(ni) = norm(cf,inf);
%     y1(ni) = norm(abscg,inf);
%     y2(ni) = norm(cl,inf);
% end
% figure(2);
% subplot(3,1,1);
% showrate(x,y0)
% subplot(3,1,2);
% showrate(x,y1)
% subplot(3,1,3);
% showrate(x,y2)


%% Test for singularity treatment right hand side
% RHS \sim \omega^0.5, the first term dominates

% h = 1/900;  a = 1/2;
% xs = 0;  ys = 0;
% [node,elem] = squaremesh([-a,a,-a,a],h); p = node;
% r = sqrt((p(:,1)-xs).^2 + (p(:,2)-ys).^2);  
% x = 2*[80 160 320 640]*pi;
% y0 = x;  y1 = x; y2 = x; 
% epsilon = 0.137;
% a = epsilon;  b = 2*epsilon;
% for ni = 1:length(x)
%     ni
%     omega = x(ni);
%     ub = 1i/4*besselh(0,1,omega*r);
%     cf_lap = cutoff_laplacian(a,b,p,xs,ys);
%     rhs2 = ub.*cf_lap;
%     rhs2(r<=a) = 0; rhs2(r>=b) = 0;
%     
%     rhs = sing_rhs(epsilon,omega,p,xs,ys);
%     rhs1 = rhs - rhs2;
%     rhs1(r<=a) = 0; rhs1(r>=b) = 0;
%     
%     y0(ni) = norm(rhs)*h;
%     y1(ni) = norm(rhs1)*h;
%     y2(ni) = norm(rhs2)*h;
% end
% figure(3);
% subplot(3,1,1);
% showrate(x,y0)
% subplot(3,1,2);
% showrate(x,y1)
% subplot(3,1,3);
% showrate(x,y2)



%% Test for the error of solution
% error \sim \omega^-0.5
% relative l2 error \sim constant

epsilon = 0.103;
x = [40 60 80 120 160 240 320]*pi;
y1 = x;
y2 = x;
y3 = x;
y4 = x;
for ni = 1:length(y1)
    ni
    omega = x(ni);
    error = error_homogeneous(omega,epsilon);
    y1(ni) = error(1);
    y2(ni) = error(2);
    y3(ni) = error(3);
    y4(ni) = error(4);
end

%%
figure(1);
subplot(2,2,1);
FJ_showrate(x,y1);
subplot(2,2,2);
FJ_showrate(x,y2);
subplot(2,2,3);
FJ_showrate(x,y3);
subplot(2,2,4);
FJ_showrate(x,y4);



