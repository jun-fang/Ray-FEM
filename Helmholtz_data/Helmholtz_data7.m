function pde = Helmholtz_data7
%% Varying slowness problem:            Need to modify the phase things!!! as well as Phase-FEM
%  uexact = exp(1i*omega/2*((x+2)^2 - (y+2)^2))

fprintf(['-'*ones(1,80) '\n'])
fprintf('Varying slowness problem: \n\n  u_ex = exp(1i*omega/2*((x+2)^2 - (y+2)^2)) \n\n');


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

    % right hand side function
    function rhs =  f(p)
        rhs = 0*p(:,1);
    end

    % speed
    function c = speed(p)
        x = p(:,1);
        y = p(:,2);
        c = 1./sqrt((x+2).^2 + (y+2).^2);   % [sqrt(2)/3, sqrt(2)/5]
    end

    % exact solution
    function u =  ex_u(p)
        global omega
        x = p(:,1); y = p(:,2);
        u = exp(1i*omega/2*((x+2).^2 - (y+2).^2));
    end
   
    % number of phases
    pde.Nphase = 1;  
 
    % exact ray
    function ray =  ray(p)
        xx = p(:,1)+2;  yy = -p(:,2)-2;
        rr = sqrt(xx.^2 + yy.^2);
        ray_dir = [xx./rr, yy./rr];
        exang = atan2(ray_dir(:,2), ray_dir(:,1));
        ray = exp(1i*exang);
    end

    % ray direction angles
    function ray_ang = ray_ang(p)
        xx = p(:,1)+2;  yy = -p(:,2)-2;
        rr = sqrt(xx.^2 + yy.^2);
        ray_dir = [xx./rr, yy./rr];
        exang = atan2(ray_dir(:,2), ray_dir(:,1));
        ray_ang = exang + 2*pi*(exang<0);
    end
    
    % exact phase
    function phase =  phase(p)
        x = p(:,1); y = p(:,2);
        phase = 1/2*((x+2).^2 - (y+2).^2);
    end
    
    % exact gradient phase
    function gradphase =  gradphase(p)
        x = p(:,1); y = p(:,2);
        gradphase = [(x+2), - (y+2)];
    end

    % gradient u
    function Du = Du(p)
        global omega
        x = p(:,1); y = p(:,2);
        exp_phase = exp(1i*omega/2*( (x+2).^2 - (y+2).^2 ));
        
        Dux = 1i*omega*(x+2).*exp_phase;
        Duy = -1i*omega*(y+2).*exp_phase;
        Du = [Dux,Duy];
    end

    % Neumann Boundary Du \dot n
    function uN =  Du_n(p)
        global a
        x = p(:,1); y = p(:,2);
        Duu = Du(p);
        Dux = Duu(:,1);
        Duy = Duu(:,2);
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

    % Impedance Boundary condition
    function g = g_IBC(p)
        global omega
        g = Du_n(p) + 1i*omega./speed(p).*ex_u(p);
    end

end