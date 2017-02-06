function pde = Helmholtz_data_four_point_sources
%% PDE data for four point sources Helmholtz problem in 2D such that 
%  
%   u_exact = sqrt(k)*besselh(0,1,k*sqrt((x+xs)^2 + (y+ys)^2))
%           + sqrt(k)*2*besselh(0,1,k*sqrt((x-xs)^2 + (y-ys)^2))
%           + sqrt(k)*0.5*besselh(0,1,k*sqrt((x+xs)^2 + (y-ys)^2))
%           - sqrt(k)*besselh(0,1,k*sqrt((x-xs)^2 + (y+ys)^2))


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
        x = p(:,1);   y = p(:,2);  k = omega;
        u = besselh(0,1,k*sqrt((x+xs).^2 + (y+ys).^2))...
            + 2*besselh(0,1,k*sqrt((x-xs).^2 + (y-ys).^2))...
            + 0.5*besselh(0,1,k*sqrt((x+xs).^2 + (y-ys).^2))...
            - besselh(0,1,k*sqrt((x-xs).^2 + (y+ys).^2));
        u = sqrt(k)*u;
    end

% exact ray
    function ray = ray(p,xs,ys)
        x = p(:,1);  y = p(:,2);
        ray_ang = zeros(length(x),4);
        ray_ang(:,1) = atan2(y+ys,x+xs);
        ray_ang(:,2) = atan2(y+ys,x-xs);
        ray_ang(:,3) = atan2(y-ys,x-xs);
        ray_ang(:,4) = atan2(y-ys,x+xs);
        ray = exp(1i*ray_ang);
    end

% gradient u
    function Du =  Du(p,xs,ys, omega)
        k = omega;
        x = p(:,1); y = p(:,2);
        r1 = sqrt((x+xs).^2 + (y+ys).^2);
        r2 = sqrt((x-xs).^2 + (y+ys).^2);
        r3 = sqrt((x-xs).^2 + (y-ys).^2);
        r4 = sqrt((x+xs).^2 + (y-ys).^2);
        
        Dux = - besselh(1,1,k*r1).*k.*(x+xs)./r1...
            + besselh(1,1,k*r2).*k.*(x-xs)./r2...
            - 2*besselh(1,1,k*r3).*k.*(x-xs)./r3...
            - 0.5*besselh(1,1,k*r4).*k.*(x+xs)./r4;
        Duy = - besselh(1,1,k*r1).*k.*(y+ys)./r1...
            + besselh(1,1,k*r2).*k.*(y+ys)./r2...
            - 2*besselh(1,1,k*r3).*k.*(y-ys)./r3...
            - 0.5*besselh(1,1,k*r4).*k.*(y-ys)./r4;
        
        Du = [sqrt(k)*Dux,sqrt(k)*Duy];
    end

% Neumann Boundary Du \dot n,  assume the domain is [-a,a]^2
    function uN =  Du_n(p,xs,ys,omega, a) 
        x = p(:,1); y = p(:,2);
        DDu = Du(p,xs,ys,omega);
        Dux = DDu(:,1);   Duy = DDu(:,2);
        uN = Dux;
        
        east_idx = find(abs(x-a)<eps);
        uN(east_idx) = Dux(east_idx);
        west_idx = find(abs(x+a)<eps);
        uN(west_idx) = -Dux(west_idx);
        
        north_idx = find(abs(y-a)<eps);
        uN(north_idx) = Duy(north_idx);
        south_idx = find(abs(y+a)<eps);
        uN(south_idx) = -Duy(south_idx);
    end

% Impedance Boundary Condition
    function g =  g_IBC(p,xs,ys, omega, a)
        g = Du_n(p,xs,ys,omega, a) + 1i*omega*u_ex(p,xs,ys, omega);
    end
end