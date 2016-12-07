function s = basic_cutoff(t)
%% basic smooth cut-off function 
pow = 2*exp(-1./abs(t))./(t-1);
s = exp(pow);

s(t<=0) = 1; s(t>=1) = 0;