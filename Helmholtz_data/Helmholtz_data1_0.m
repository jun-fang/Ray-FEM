function pde = Helmholtz_data1_0
%% uexact = sqrt(k)*first Hankel(0,1,k*sqrt((x-xs)^2 + (y-ys)^2))

fprintf(['-'*ones(1,80) '\n'])
fprintf('One point source problem (Outside domain): \n\n  u_ex = sqrt(k)*besselh(0,1,k*sqrt((x-2)^2 + (y-2)^2)) \n\n');

pde = struct(...
    'f',@f,...                    % right hand side source
    'speed',@speed,...            % medium speed
    'ex_u',@ex_u,...              % exact solution
    'Nphase',@Nphase,...          % number of phases
    'phase',@phase,...            % exact phase
    'gradphase',@gradphase,...    % gradient of exact phase 
    'int_coe',@int_coe,...        % nodal interpolation coefficient
    'ray',@ray,...                % exact ray direction
    'Du',@Du,...                  % gradient of exact solution  
    'Du_n',@Du_n,...              % Neumann boundary condition \partial u / \partial n or Du \cdot n
    'g_IBC',@g_IBC);              % Impedance boundary condition  \partial u / \partial n + iku = g_IBC


global omega;
global xs;
global ys;


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
        k = omega;
        x = p(:,1); y = p(:,2);
        u = sqrt(k)*besselh(0,1,k*sqrt((x-xs).^2 + (y-ys).^2));
    end
        
    % number of phases
    pde.Nphase = 1;
    
    % exact phase
    function phase =  phase(p)
        xx = p(:,1)-xs; yy = p(:,2)-ys;
        phase = sqrt(xx.^2 + yy.^2);
    end

    % gradient phase
    function gphase =  gradphase(p)
        xx = p(:,1)-xs; yy = p(:,2)-ys;
        ray_dir = atan2(yy,xx);
        gphase = exp(1i*ray_dir);  
    end

    % nodal interpolation coefficient 
    function int_coe = int_coe(p)
        k = omega;
        x = p(:,1);   y = p(:,2);
        amp = ex_u(p)./exp(1i*k*phase(p));
        gpha = gradphase(p);
        temp = real(gpha).*x + imag(gpha).*y;
        pha = phase(p) - temp;
        int_coe = amp.*exp(1i*k*pha); 
    end

    % exact ray
    function ray = ray(p)
        xx = p(:,1)-xs; yy = p(:,2)-ys;
        ray = atan2(yy,xx);
        ray = ray + 2*pi*(ray<0);
        
%         ray = exp(1i*ray);
    end

    % gradient u
    function Du =  Du(p)
        k = omega;
        x = p(:,1)-xs; y = p(:,2)-ys;
        r = sqrt(x.^2 + y.^2);
        
        Dux = -sqrt(k)*besselh(1,1,k*r).*k.*x./r;
        Duy = -sqrt(k)*besselh(1,1,k*r).*k.*y./r;
        Du = [Dux, Duy];
    end

    % Neumann Boundary Du \dot n
    function uN =  Du_n(p)
        global a
        k = omega;
        x = p(:,1)-xs; y = p(:,2)-ys;
        r = sqrt(x.^2 + y.^2);
        xx = p(:,1);  yy = p(:,2);
        
        Dux = -sqrt(k)*besselh(1,1,k*r).*k.*x./r;
        Duy = -sqrt(k)*besselh(1,1,k*r).*k.*y./r;
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
        k = omega;
        g = Du_n(p) + 1i*k*ex_u(p);
    end
end