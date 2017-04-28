function f = arccosh(z, der)
%% Inverse of hyperpolic cosine function

if der == 0 
    f = log(z + sqrt(z.*z - 1));
end

if der == 1   % first derivative
    f = 1./sqrt(z.*z - 1);
end
