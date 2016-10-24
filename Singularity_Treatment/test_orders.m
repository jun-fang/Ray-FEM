%% Test for asymptotic order of hankel function 
% H(z) \sim O(z^-0.5)

% n = [1,2,3,4,5];
% x = 2.^n;
% y1 = besselh(0,1,x);
% y2 = besselh(1,1,x);
% showrate(x,abs(y1))
% showrate(x,abs(y2))


%% Test for cut-off function wrt epsilon
% gradient \sim eps^0.5,
% laplacian \sim eps^2.5

% h = 1/3900;  a = 1/2;
% [node,elem] = squaremesh([-a,a,-a,a],h);
% p = node;
% n = [1,2,3,4,5];
% x = 1.7.^(n+2);
% y0 = x; y1 = x; y2 = x;
% for ni = 1:length(n)
%     ni
%     epsilon = x(ni);
%     a = 1/epsilon;  b = 2/epsilon;
%     cg = cutoff_gradient(a,b,p);
% %     abscg = sqrt(abs(cg(:,1)).^2 + abs(cg(:,1)).^2);
% %     cl = cutoff_laplacian(a,b,p);
%     cf = cutoff(a,b,p);
%     y0(ni) = norm(cf);
% %     y1(ni) = norm(abscg);
% %     y2(ni) = norm(cl);
% end
% showrate(x,y0);
% showrate(x,y1)
% showrate(x,y2)


%% Test for singularity treatment right hand side
% RHS \sim \omega^0.5

% h = 1/900;  a = 1/2;
% [node,elem] = squaremesh([-a,a,-a,a],h); p = node;
% r = sqrt(p(:,1).^2 + p(:,2).^2);  
% x = 2*[80 160 320 640]*pi;
% y = x;  y1 = x; y2 = x; 
% epsilon = 1/0.13;
% a = 1/epsilon;  b = 2/epsilon;
% for ni = 1:length(y)
%     ni
%     omega = x(ni);
% %     f = sing_rhs(epsilon,omega,p);
% %     y(ni) = norm(f)*h;
%     
%     f1 = -1i/4*omega*besselh(1,1,omega*r)./(r.*r).*p(:,1);
%     f1(r<a) = 0;  f1(r>b) = 0;
%     y1(ni) = norm(f1)*h;
% end
% showrate(x,y1)


%% Test for the error of solution
% error \sim \omega^-0.5
% relative l2 error \sim constant

epsilon = 0.137;
x = [40 60 80 120 160 240 320 480 640]*pi;
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
figure(1);
subplot(2,2,1);
FJ_showrate(x,y1);
subplot(2,2,2);
FJ_showrate(x,y2);
subplot(2,2,3);
FJ_showrate(x,y3);
subplot(2,2,4);
FJ_showrate(x,y4);



