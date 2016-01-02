function pde = Helmholtz_data4
%% uexact = sqrt(k)*besselh(0,1,k*sqrt((x+2)^2 + (y+2)^2))
%         + exp(1i*k*sqrt(2)/2*(x-y))

fprintf(['-'*ones(1,80) '\n'])
fprintf('Two phases problem: \n\n  u_ex = sqrt(k)*besselh(0,1,k*sqrt((x+2)^2 + (y+2)^2)) + exp(1i*k*sqrt(2)/2*(x-y)) \n\n');

pde = struct(...
    'f',@f,...                    % right hand side source
    'speed',@speed,...            % medium speed
    'ex_u',@ex_u,...              % exact solution
    'Nphase',@Nphase,...          % number of phases
    'phase',@phase,...            % exact phase
    'gradphase',@gradphase,...    % gradient of exact phase 
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
        global omega;
        k = omega;
        x = p(:,1);    y = p(:,2);
        u = sqrt(k)*besselh(0,1,k*sqrt((x+2).^2 + (y+2).^2)) + exp(1i*k*sqrt(2)/2*(x-y));
    end

    % number of phases
    pde.Nphase = 2;

    % exact phase
    function phase =  phase(p)
        N = length(p(:,1));
        phase = zeros(N,2);
        x = p(:,1);  y = p(:,2);
        
        phase(:,1) = sqrt((x+2).^2 + (y+2).^2);
        phase(:,2) = sqrt(2)/2*(x-y);
    end

    % gradient phase
    function gphase =  gradphase(p)
        N = length(p(:,1));
        angle = zeros(N,2);
        x = p(:,1);  y = p(:,2);
        
        angle(:,1) = atan2(y+2,x+2);
        angle(:,2) = atan2(-ones(N,1), ones(N,1));
        gphase = exp(1i*angle);  
    end

    % exact ray
    function ray = ray(p)
        ray = gradphase(p);
    end

    % gradient u
    function Du =  Du(p)
        global omega;
        k = omega;
        x = p(:,1); y = p(:,2);
        r = sqrt((x+2).^2 + (y+2).^2);
        
        Du2 = 1i*k*sqrt(2)/2*exp(1i*k*sqrt(2)/2*(x-y));
        
        Dux = -sqrt(k)*besselh(1,1,k*r).*k.*(x+2)./r + Du2;
        Duy = -sqrt(k)*besselh(1,1,k*r).*k.*(y+2)./r - Du2;
        Du = [Dux,Duy];
    end

    % Neumann Boundary Du \dot n  
    function uN =  Du_n(p)
        global omega;
        global a
        k = omega;
        x = p(:,1); y = p(:,2);
        r = sqrt((x+2).^2 + (y+2).^2);
        
        Du2 = 1i*k*sqrt(2)/2*exp(1i*k*sqrt(2)/2*(x-y));
        
        Dux = -sqrt(k)*besselh(1,1,k*r).*k.*(x+2)./r + Du2;
        Duy = -sqrt(k)*besselh(1,1,k*r).*k.*(y+2)./r - Du2;
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
        global omega;
        k = omega;
        g = Du_n(p) + 1i*k*ex_u(p);
    end
end