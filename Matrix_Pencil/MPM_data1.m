function pde = MPM_data1
%% u_exact = sqrt(omega)*first Hankel(0,1,omega*sqrt((x-xs)^2 + (y-ys)^2))

global omega xs ys a;

pde = struct(...
    'f',@f,...                    % right hand side source
    'speed',@speed,...            % medium speed
    'ex_u',@ex_u,...              % exact solution
    'ray',@ray,...                % exact ray direction
    'Du',@Du,...                  % gradient of exact solution  
    'Du_n',@Du_n,...              % Neumann boundary condition \partial u / \partial n or Du \cdot n
    'g_IBC',@g_IBC);              % Impedance boundary condition  \partial u / \partial n + iku = g_IBC

    % right hand side 
    function rhs =  f(p)
        rhs = 0*p(:,1);
    end

    % medium speed
    function c = speed(p)
        n = length(p(:,1));
        c = ones(n,1);
    end

    % exact solution
    function u =  ex_u(p)
        x = p(:,1); y = p(:,2);
        u = sqrt(omega)*besselh(0,1,omega*sqrt((x-xs).^2 + (y-ys).^2));
    end
        
    % exact ray
    function ray = ray(p)
        xx = p(:,1)-xs; yy = p(:,2)-ys;
        ray = atan2(yy,xx);
        ray = ray + 2*pi*(ray<0);
        rr = sqrt(xx.^2 + yy.^2);
        ray = ray.*(rr>10*eps);
    end

    % gradient u
    function Du =  Du(p)
        x = p(:,1)-xs; y = p(:,2)-ys;
        r = sqrt(x.^2 + y.^2);
        
        Dux = -sqrt(omega)*besselh(1,1,omega*r).*omega.*x./r;
        Duy = -sqrt(omega)*besselh(1,1,omega*r).*omega.*y./r;
        Du = [Dux, Duy];
    end

    % Neumann Boundary Du \dot n
    function uN =  Du_n(p)
        DDu = Du(p);  uN = DDu(:,1);
        Dux = DDu(:,1);  Duy = DDu(:,2);
        xx = p(:,1);  yy = p(:,2);
        
        east_idx = find(abs(xx-a)<eps);
        uN(east_idx) = Dux(east_idx);
        west_idx = find(abs(xx+a)<eps);
        uN(west_idx) = -Dux(west_idx);
        
        north_idx = find(abs(yy-a)<eps);
        uN(north_idx) = Duy(north_idx);
        south_idx = find(abs(yy+a)<eps);
        uN(south_idx) = -Duy(south_idx);
    end

    % Impedance Boundary Condition
    function g =  g_IBC(p)
        g = Du_n(p) + 1i*omega*ex_u(p);
    end

end