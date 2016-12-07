function s1 = basic_cutoff_first_derivative(t)
%% The first derivative of the basic smooth cut-off function 
pow = 2*exp(-1./abs(t))./(t-1);
pow_der = 2*exp(-1./abs(t)).*(t-1-t.^2)...
    ./( t.*(t-1) ).^2;
s1 = exp(pow).*pow_der;

s1(t<=0) = 0; s1(t>=1) = 0;