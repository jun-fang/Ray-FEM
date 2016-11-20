function pde = Helmholtz_data_point_source
%% PDE data for Green function of Helmholtz Equation in 2D
% PED:      - (\Delta + \omega^2)u = \delta([x-xs, y-ys])
% Solution: uexact = i/4*first Hankel(0,1,omega*sqrt((x-xs)^2 + (y-ys)^2))

global omega xs ys;


fprintf(['-'*ones(1,80) '\n'])
fprintf('Point source problem:    (xs,ys) = (%d,%d)\n\n  ',xs,ys);
fprintf('u_ex = 1i/4*besselh(0,1,omega*sqrt((x-xs)^2 + (y-ys)^2)) \n\n');


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



% right hand side
    function rhs =  f(p)
        x = p(:,1);  y = p(:,2);
        r = sqrt((x-xs).^2 + (y-ys).^2);
        rhs = r<eps;
    end

% medium speed
    function c = speed(p)
        c = 0*p(:,1)+1;
    end

% exact solution
    function u =  ex_u(p)
        x = p(:,1); y = p(:,2);
        u = 1i/4*besselh(0,1,omega*sqrt((x-xs).^2 + (y-ys).^2));
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
        gphase = ray(p);
    end

% nodal interpolation coefficient
    function int_coe = int_coe(p)
        x = p(:,1);   y = p(:,2);
        amp = ex_u(p)./exp(1i*omega*phase(p));
        gpha = gradphase(p);
        temp = real(gpha).*x + imag(gpha).*y;
        pha = phase(p) - temp;
        int_coe = amp.*exp(1i*omega*pha);
    end

% exact ray
    function ray = ray(p)
        xx = p(:,1)-xs; yy = p(:,2)-ys;
        ray = atan2(yy,xx);
        ray = exp(1i*ray);
    end

% gradient u
    function Du =  Du(p)
        x = p(:,1)-2; y = p(:,2)-2;
        r = sqrt(x.^2 + y.^2);
        
        Dux = -1i/4*besselh(1,1,omega*r).*omega.*x./r;
        Duy = -1i/4*besselh(1,1,omega*r).*omega.*y./r;
        Du = [Dux, Duy];
    end

% Neumann Boundary Du \dot n
    function uN =  Du_n(p)
        global a;
        x = p(:,1)-2; y = p(:,2)-2;
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
    function g =  g_IBC(p)
        g = Du_n(p) + 1i*omega*ex_u(p);
    end
end