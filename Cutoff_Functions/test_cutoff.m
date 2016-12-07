
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


%% test the order of cut-off, its gradient and Laplacian with respect to epsilon

