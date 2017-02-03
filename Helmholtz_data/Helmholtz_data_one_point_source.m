function pde = Helmholtz_data_one_point_source
%% PDE data for point source Helmholtz problem in 2D
%  
%         PDE:      - (\Delta + \omega^2)u = \delta([x-xs, y-ys])
%
%   Exact solution: u_exact = 1i/4*besselh(0,1,omega*sqrt((x-xs)^2 + (y-ys)^2))
%

pde = struct(...
    'f',@f,...                    % right hand side source
    'speed',@speed,...            % medium speed
    'u_ex',@u_ex,...              % exact solution
    'ray',@ray,...                % exact ray direction
    'Du',@Du,...                  % gradient of exact solution
    'Du_n',@Du_n,...              % Neumann boundary condition \partial u / \partial n or Du \cdot n
    'g_IBC',@g_IBC);              % Impedance boundary condition  \partial u / \partial n + iku = g_IBC



% right hand side
    function rhs =  f(p,xs,ys)
        x = p(:,1);  y = p(:,2);
        r = sqrt((x-xs).^2 + (y-ys).^2);
        rhs = r<eps;
    end

% medium speed
    function c = speed(p)
        c = 0*p(:,1)+1;
    end

% exact solution
    function u =  u_ex(p,xs,ys, omega)
        x = p(:,1); y = p(:,2);
        u = 1i/4*besselh(0,1,omega*sqrt((x-xs).^2 + (y-ys).^2));
    end

% exact ray
    function ray = ray(p,xs,ys)
        xx = p(:,1)-xs; yy = p(:,2)-ys;
        rr = sqrt(xx.*xx + yy.*yy);
        ray = atan2(yy,xx);
        ray = exp(1i*ray).*(rr>eps);
    end

% gradient u
    function Du =  Du(p,xs,ys, omega)
        x = p(:,1)-xs; y = p(:,2)-ys;
        r = sqrt(x.^2 + y.^2);
        
        Dux = -1i/4*besselh(1,1,omega*r).*omega.*x./r;
        Duy = -1i/4*besselh(1,1,omega*r).*omega.*y./r;
        Du = [Dux, Duy];
    end

% Neumann Boundary Du \dot n,  assume the domain is [-a,a]^2
    function uN =  Du_n(p,xs,ys, omega, a) 
        x = p(:,1)-xs; y = p(:,2)-ys;
        r = sqrt(x.^2 + y.^2);
        xx = p(:,1);  yy = p(:,2);
        
        Dux = -1i/4*besselh(1,1,omega*r).*omega.*x./r;
        Duy = -1i/4*besselh(1,1,omega*r).*omega.*y./r;
        uN = Dux;
        
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
    function g =  g_IBC(p,xs,ys, omega, a)
        g = Du_n(p,xs,ys, omega, a) + 1i*omega*u_ex(p,xs,ys, omega);
    end
end