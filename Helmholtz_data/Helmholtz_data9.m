function pde = Helmholtz_data9

% uexact = sqrt(k)*first Hankel(0,1,k*sqrt((x-1/8)^2 + (y-1/10)^2))

fprintf(['-'*ones(1,80) '\n'])
fprintf('One point source problem (inside domain): \n\n  u_ex = sqrt(k)*besselh(0,1,k*sqrt((x-1/8)^2 + (y-1/10)^2)) \n\n');

global omega a;



pde = struct(...
    'f',@f,...                    % right hand side source
    'speed',@speed,...            % medium speed
    'ex_u',@ex_u,...              % exact solution
    'Nphase',@Nphase,...          % number of phases
    'phase',@phase,...            % exact phase
    'gradphase',@gradphase,...    % gradient of exact phase 
    'ray',@ray,...                % exact ray direction
    'ray_ang',@ray_ang,...        % exact ray direction angles [0,2pi]
    'Du',@Du,...                  % gradient of exact solution  
    'Du_n',@Du_n,...              % Neumann boundary condition \partial u / \partial n or Du \cdot n
    'g_IBC',@g_IBC);              % Impedance boundary condition  \partial u / \partial n + iku = g_IBC

    % right hand side 
    function rhs =  f(p)
        rhs = -(1+1i*sqrt(2)*omega)*ex_u(p);
    end

    % medium speed
    function c = speed(p)
        n = length(p(:,1));
        c = ones(n,1);
    end

    % exact solution
    function u =  ex_u(p)
        u =  ( exp(p(:,1))+exp(p(:,2)) ).*exp( 1i*omega*sqrt(2)/2*( p(:,1)+p(:,2) ) );
    end
        
    % number of phases
    pde.Nphase = 1;
    
    % exact phase
    function phase =  phase(p)
        phase = sqrt(2)/2*( p(:,1)+p(:,2) );
    end

    % gradient phase
%     function gphase =  gradphase(p)
%         xx = p(:,1)-xc; yy = p(:,2)-yc;
%         ray_dir = atan2(yy,xx);
%         gphase = exp(1i*ray_dir);  
%     end

    % exact ray
%     function ray = ray(p)
%         xx = p(:,1)-xc; yy = p(:,2)-yc;
%         ray = atan2(yy,xx);
%         ray = exp(1i*ray);
%         ray = ray.*(1 - (abs(xx)<eps).*(abs(yy)<eps));
%     end
% 
%     % exact ray angle
%     function ray_ang = ray_ang(p)
%         xx = p(:,1)-xc; yy = p(:,2)-yc;
%         ray = atan2(yy,xx);
%         ray_ang = ray + 2*pi*(ray<0);
%         ray = ray.*(1 - (abs(xx)<eps).*(abs(yy)<eps));
%     end

    % gradient u
    function Du =  Du(p)
        x = p(:,1);  y = p(:,2);
        exp_phase = exp( 1i*omega*sqrt(2)/2*( x+y ) );
        temp = 1i*omega*sqrt(2)/2*( exp(x)+exp(y) );
        
        Dux = (exp(x) + temp).*exp_phase;
        Duy = (exp(y) + temp).*exp_phase;
        Du = [Dux, Duy];
    end

    % Neumann Boundary Du \dot n
    function uN =  Du_n(p)
        xx = p(:,1);  yy = p(:,2); x = xx; y = yy;
        exp_phase = exp( 1i*omega*sqrt(2)/2*( x+y ) );
        temp = 1i*omega*sqrt(2)/2*( exp(x)+exp(y) );
        
        Dux = (exp(x) + temp).*exp_phase;
        Duy = (exp(y) + temp).*exp_phase;
        Du = [Dux, Duy];
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
    function g =  g_IBC(p)
        g = Du_n(p) + 1i*omega*ex_u(p);
    end
end