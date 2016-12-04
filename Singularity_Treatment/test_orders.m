%% Test for asymptotic order of hankel function 
% H(z) \sim O(z^-0.5)

% n = [1,2,3,4,5];
% x = 2.^n;
% y1 = besselh(0,1,x);
% y2 = besselh(1,1,x);
% showrate(x,abs(y1))
% showrate(x,abs(y2))


%% Test for cut-off function wrt epsilon
% cut-off L2: \sim eps^1, L^{\infy} \sim eps^0
% gradient L2: \sim eps^-0.5, L^{\infty} \sim eps^-2
% laplacian L2: \sim eps^-2.5, L^{\infty} \sim eps^-4

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
%     abscg = sqrt(abs(cg(:,1)).^2 + abs(cg(:,1)).^2);
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
% showrate(x,y0)
% subplot(3,1,2);
% showrate(x,y1)
% subplot(3,1,3);
% showrate(x,y2)
% 

%% Test for singularity treatment right hand side
% fixed epsilon, RHS \sim \omega^0.5, the first term dominates
% with fixed omega, why the order of the second term is different from
% cutoff laplacian?? I see, the support of cut-off function depends on
% epsilon.   
%     RHS1: L^2 \sim eps^-1, L^{\infty} \sim eps^-2.5
%     RHS2: L^2 \sim eps^-3, L^{\infty} \sim eps^-4.5
%     RHS:  L^2 \sim eps^-2.2, L^{\infty} \sim eps^-3.6
% need to figure it out

h = 1/900;  a = 1/2;
xs = 0;  ys = 0;
[node,elem] = squaremesh([-a,a,-a,a],h); p = node;
r = sqrt((p(:,1)-xs).^2 + (p(:,2)-ys).^2);  
x = 1/4*[80 160 320 640]*pi;

% omega = 640*pi;
% n = [1,2,3,4,5];
% x = 1.7.^-(n+2);

y0 = x;  y1 = x; y2 = x; 
epsilon = 0.137;
a = epsilon;  b = 2*epsilon;
for ni = 1:length(x)
    ni
%     epsilon = x(ni);
%     a = epsilon;  b = 2*epsilon;

    omega = x(ni);
    ub = 1i/4*besselh(0,1,omega*r);
    cf_lap = cutoff_laplacian(a,b,p,xs,ys);
    rhs2 = ub.*cf_lap;
    rhs2(r<=a) = 0; rhs2(r>=b) = 0;
    
    rhs = sing_rhs_homo(epsilon,omega,p,xs,ys);
    rhs(r<=a) = 0; rhs(r>=b) = 0;
    rhs1 = rhs - rhs2;
    rhs1(r<=a) = 0; rhs1(r>=b) = 0;
    
    y0(ni) = norm(rhs,inf)*h;
    y1(ni) = norm(rhs1,inf)*h;
    y2(ni) = norm(rhs2,inf)*h;
end
figure(3);
subplot(3,1,1);
showrate(x,y0)
subplot(3,1,2);
showrate(x,y1)
subplot(3,1,3);
showrate(x,y2)



%% Test for the error of solution
% error \sim \omega^-0.5
% relative l2 error \sim constant

% epsilon = 0.103;
% x = pi*[40 60 80 120 160]% 240 320];
% epss = 1/2*x.^(-.2)
% y1 = x;
% y2 = x;
% y3 = x;
% y4 = x;
% y0 = x;
% for ni = 1:length(y1)
%     ni
%     omega = x(ni);
%     epsilon = epss(ni);
%     error = error_homogeneous(omega,epsilon);
%     
%     y1(ni) = error(1);
%     y2(ni) = error(2);
%     y3(ni) = error(3);
%     y4(ni) = error(4);
%     y0(ni) = error(5);
% end

%%
% figure(1);
% subplot(2,2,1);
% FJ_showrate(x,y1);
% subplot(2,2,2);
% FJ_showrate(x,y2);
% subplot(2,2,3);
% FJ_showrate(x,y3);
% subplot(2,2,4);
% FJ_showrate(x,y4);
% 
% figure(2);
% FJ_showrate(x,y0);
% 


