
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

h = 1/1900;  a = 1/2;
[node,elem] = squaremesh([-a,a,-a,a],h);
p = node;
n = [1,2,3,4,5];
x = 1.7.^-(n+2);
y0 = x; y1 = x; y2 = x;
for ni = 1:length(n)
    ni
    epsilon = x(ni);
    a = epsilon;  b = 2*epsilon;
    cg = cutoff_gradient(p,0,0,a,b);
    abscg = abs(cg(:,1));%sqrt(abs(cg(:,1)).^2 + abs(cg(:,1)).^2);
    cl = cutoff_laplacian(p,0,0,a,b);
    cf = cutoff(p,0,0,a,b);
    y0(ni) = norm(cf);
    y1(ni) = norm(abscg);
    y2(ni) = norm(cl);
%     y0(ni) = norm(cf,inf);
%     y1(ni) = norm(abscg,inf);
%     y2(ni) = norm(cl,inf);
end
figure(2);
subplot(3,1,1);
show_convergence_rate(x,y0)
subplot(3,1,2);
show_convergence_rate(x,y1)
subplot(3,1,3);
show_convergence_rate(x,y2)
