function pde = Helmholtz_data2_0
%% uexact = exp(1i*k*sqrt(2)/2*(x+y))

fprintf(['-'*ones(1,80) '\n'])
fprintf('Perfect plane wave problem: \n\n  u_ex = exp(1i*k*sqrt(2)/2*(x+y)) \n\n');


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


    

    % medium speed
    function c = speed(p)
        n = length(p(:,1));
        c = ones(n,1);
    end 

    % exact solution
    function u =  ex_u(p)
        global omega;
        k = omega;
        x = p(:,1); y = p(:,2);
        u = exp(1i*k*sqrt(2)/2*(x+y));
    end

    % number of phases
    pde.Nphase = 1;
    
    
    % exact phase
    function phase =  phase(p)
        x = p(:,1); y = p(:,2); 
        phase = sqrt(2)/2*(x+y);
    end

    % gradient phase
    function gphase =  gradphase(p)
        n = length(p(:,1));
        gphase = exp(1i*pi/4)*ones(n,1);  
    end

    % exact ray
    function ray = ray(p)
        n = length(p(:,1));
        ray = exp(1i*pi/4)*ones(n,1);  
    end

    % gradient u
    function Du =  Du(p)
        global omega;
        k = omega;
        x = p(:,1); y = p(:,2);
        
        exp_phase = exp(1i*k*sqrt(2)/2*(x+y));
        Ux = 1i*k*sqrt(2)/2.*exp_phase;
        Uy = 1i*k*sqrt(2)/2.*exp_phase;
        Du = [Ux, Uy];
    end
    
    % Neumann Boundary Du \dot n   
    function uN =  Du_n(p)
        global omega;
        global a;
        k = omega;
        x = p(:,1); y = p(:,2);
        
        exp_phase = exp(1i*k*sqrt(2)/2*(x+y));
        Ux = (exp(x) + 1i*k*sqrt(2)/2*(exp(x) + exp(y))).*exp_phase;
        Uy = (exp(y) + 1i*k*sqrt(2)/2*(exp(x) + exp(y))).*exp_phase;
        uN = exp_phase;
        
        east_idx = find(abs(x-a)<eps);
        uN(east_idx) = Ux(east_idx);
        
        west_idx = find(abs(x+a)<eps);
        uN(west_idx) = -Ux(west_idx);
        
        north_idx = find(abs(y-a)<eps);
        uN(north_idx) = Uy(north_idx);
        
        south_idx = find(abs(y+a)<eps);
        uN(south_idx) = -Uy(south_idx);
    end

    % Impedance Boundary Condition
    function g =  g_IBC(p)
        global omega;
        k = omega;
        g = Du_n(p) + 1i*k*ex_u(p);
    end    
    
end