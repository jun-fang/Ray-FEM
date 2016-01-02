function pde = Helmholtz_data6
%% Gaussina distributed source problem
% source term f = exp(-100*((x-1/4)^2 + (y-1/8)^2))

fprintf(['-'*ones(1,80) '\n'])
fprintf('Gaussina distributed source Problem: \n\n  f = exp(-100*((x-1/4)^2 + (y-1/8)^2)) \n\n');


pde = struct('f',@f,'speed',@speed);

    % load data (right hand side function)
    function rhs =  f(p)
        x = p(:,1);    y = p(:,2);
        rhs = exp(-100*((x-1/4).^2 + (y-1/8).^2));
    end

    function speed = speed(p)
        n = length(p(:,1));
        speed = ones(n,1);
    end

end