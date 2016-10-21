function pde = MPM_data2_0
%% uexact = sqrt(k)*besselh(0,1,k*sqrt((x+xs)^2 + (y+ys)^2))
%         - sqrt(k)*besselh(0,1,k*sqrt((x-xs)^2 + (y-ys)^2))


global omega xs ys a;

pde = struct(...
    'f',@f,...                    % right hand side source
    'speed',@speed,...            % medium speed
    'ex_u',@ex_u,...              % exact solution
    'Nphase',@Nphase,...          % number of phases
    'phase',@phase,...            % exact phase
    'gradphase',@gradphase,...    % gradient of exact phase 
    'ray',@ray,...                % exact ray information, stored as complex form e^{i*\theta}
    'ray_ang',@ray_ang,...        % exact ray direction angles [0,2pi]
    'Du',@Du,...                  % gradient of exact solution  
    'Du_n',@Du_n,...              % Neumann boundary condition \partial u / \partial n or Du \cdot n
    'g_IBC',@g_IBC);              % Impedance boundary condition  \partial u / \partial n + iku = g_IBC

    % right hand side function
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
        k = omega;
        x = p(:,1);    y = p(:,2);
        u = besselh(0,1,k*sqrt((x+xs).^2 + (y+ys).^2))...
            - besselh(0,1,k*sqrt((x-xs).^2 + (y-ys).^2));
        u = sqrt(k)*u;
    end

    % number of phases
    pde.Nphase = 2;

    % exact phase
    function phase =  phase(p)
        N = length(p(:,1));
        phase = zeros(N,2);
        x = p(:,1);  y = p(:,2);
        
        phase(:,1) = sqrt((x+xs).^2 + (y+ys).^2);
        phase(:,2) = sqrt((x-xs).^2 + (y-ys).^2);
    end

    % ray direction angles
    function ray_ang = ray_ang(p)
        N = length(p(:,1));
        ray_ang = zeros(N,2);
        x = p(:,1);  y = p(:,2);
        
        ray_ang(:,1) = atan2(y+ys,x+xs);
        ray_ang(:,2) = atan2(y-ys,x-xs);
        ray_ang = ray_ang + 2*pi*(ray_ang<0);
    end

    % gradient phase
    function gphase =  gradphase(p)
        angle = ray_ang(p);
        gphase = exp(1i*angle);  
    end

    % exact ray
    function ray = ray(p)
        ray = ray_ang(p);
    end

    % gradient u
    function Du =  Du(p)
        k = omega;
        x = p(:,1); y = p(:,2);
        r1 = sqrt((x+xs).^2 + (y+ys).^2);
        r2 = sqrt((x-xs).^2 + (y-ys).^2);
        
        Dux = - besselh(1,1,k*r1).*k.*(x+xs)./r1...
            + besselh(1,1,k*r2).*k.*(x-xs)./r2;

        Duy = - besselh(1,1,k*r1).*k.*(y+ys)./r1...
            + besselh(1,1,k*r2).*k.*(y-ys)./r2;

        
        Du = [sqrt(k)*Dux,sqrt(k)*Duy];
    end

    % Neumann Boundary Du \dot n  
    function uN =  Du_n(p)
        x = p(:,1); y = p(:,2);
        DDu = Du(p);
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
    function g =  g_IBC(p)
        k = omega;
        g = Du_n(p) + 1i*k*ex_u(p);
    end
end