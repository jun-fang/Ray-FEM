function s2 = basic_cutoff_second_derivative(t)
%% The second derivative of the basic smooth cut-off function 
pow = 2*exp(-1./abs(t))./(t-1);
N = 2*exp(-1./abs(t)).*(t-1-t.^2);  % numerator of the pow's derivative
D = ( t.*(t-1) ).^2;                % denominator of the pow's derivative
pow_der = N./D;                     % pow's derivative

N_der = 2*exp(-1./abs(t)).*(t-1-2*t.^3)./(t.^2);
D_der = 2*t.*(t-1).*(2*t-1);
pow_der_der = (N_der.*D - N.*D_der)./(D.^2);  % pow's second derivative

s2 = exp(pow).*( pow_der.^2 + pow_der_der);

s2(t<=0) = 0;   s2(t>=1) = 0;