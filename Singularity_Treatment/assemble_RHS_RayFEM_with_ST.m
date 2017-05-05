function b = assemble_RHS_RayFEM_with_ST(node,elem,xs,ys,omega,epsilon,wpml,sigmaMax,ray,speed,fquadorder,option)
%% Function to assemble the right hand side with singularity treatment (ST):
%         -\Delta u - (omega/c)^2 u = 2 \nabla u_b \cdot \nabla \chi + u_b \Delta \chi,     in D
%                                                u = 0               on \partial D
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT:
%
%   node: N x 2 matrix that contains the physical position of each node
%         node(:,1) provides the x coordinate
%         node(:,2) provides the y coordinate
%
%   elem: NT x 3 matrix that contains the indices of the nodes for each
%         triangle element
%
%   (xs, ys): source location;   omega: frequency;
%
%   epsilon: cut-off support size
%
%   wpml: width of PML;  sigmaMax: max absorption
%
%   ray: ray information;  speed: wave speed
%
%   fquadorder: The order of numerical quadrature
%
%   option: different types -- homogenenous, gravity
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% OUTPUT:
%
%   b: Ndof x 1 Galerkin projection of the source
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fprintf('Assembling the right-hand side \n');

% source = @(x) -1i*10*pi*exp(1i*10*pi*sqrt(x(:,1).^2+x(:,2).^2))./sqrt(x(:,1).^2+x(:,2).^2);



%% FEM set up
N = size(node,1);        % number of grid points
NT = size(elem,1);       % number of triangle elements
Nray = size(ray,2);
Ndof = N*Nray;          % degree of freedom

k = omega./speed(node);          % wavenumber
kk = repmat(k,1,Nray);

xmax = max(node(:,1));
xmin = min(node(:,1));
ymax = max(node(:,2));
ymin = min(node(:,2));
h = (xmax - xmin)/(round(sqrt(N)) - 1);


%% PML set up
sigmaPML_x = @(x)sigmaMax*( (x-xmin-wpml).^2.*(x < xmin + wpml) + ...
    (x-(xmax-wpml)).^2.*(x > xmax - wpml))/wpml^2;
sigmaPML_y = @(y) sigmaMax*( (y-ymin-wpml).^2.*(y < ymin + wpml) ...
    + (y-(ymax-wpml)).^2.*(y > ymax - wpml))/wpml^2;
s_xy = @(x,y) ((1+1i*sigmaPML_x(x)/omega).*(1+1i*sigmaPML_y(y)/omega));    %% 1/(s1*s2)


%% Numerical Quadrature
[lambda,weight] = quadpts(fquadorder);
phi = lambda;           % linear bases
nQuad = size(lambda,1);


%% Compute geometric quantities and gradient of local basis
[~,area] = gradbasis(node,elem);
reparea = repmat(area,Nray,1);


%% Assemble right-hand side

bt = zeros(NT*Nray,3);       % the right hand side
b = zeros(Ndof,1);

%% Babich pre-processing
if  iscell(option) && strcmp(option{1}, 'Babich')
    
    if strcmp(option{2}, 'CGV')
        switch round(omega/(pi))
            case 120
                load('Babich_CGV_30.mat');
            case 160
                load('Babich_CGV_40.mat');
            case 240
                load('Babich_CGV_60.mat');
            case 320
                load('Babich_CGV_80.mat');
            case 400
                load('Babich_CGV_100.mat');
            case 500
                load('Babich_CGV_125.mat');
            case 600
                load('Babich_CGV_150.mat');
        end
    end
    
    if strcmp(option{2}, 'Homo')
        switch round(omega/(pi))
            case 100
                load('Babich_Homo_25.mat');
            case 160
                load('Babich_Homo_40.mat');
            case 240
                load('Babich_Homo_60.mat');
            case 400
                load('Babich_Homo_100.mat');
            case 600
                load('Babich_Homo_150.mat');
        end
    end
    
    a = 1/2;
    CompressRatio = round(Bh0/(h/2));
    Bh = 1/round( 1/(Bh0/CompressRatio) );
    Bx = -a: Bh : a;  By = -a: Bh : a;
    [BX0, BY0] = meshgrid(Bx0, By0);
    [BX, BY] = meshgrid(Bx, By);
    
    %% refined amplitude
    DD1 = interp2(BX0,BY0,D1,BX,BY,'spline');
    DD2 = interp2(BX0,BY0,D2,BX,BY,'spline');
    
    % gradient
    [D1x,D1y] = num_derivative(D1,Bh0,4);
    [D2x,D2y] = num_derivative(D2,Bh0,4);
    DD1x = interp2(BX0,BY0,D1x,BX,BY,'spline');
    DD1y = interp2(BX0,BY0,D1y,BX,BY,'spline');
    DD2x = interp2(BX0,BY0,D2x,BX,BY,'spline');
    DD2y = interp2(BX0,BY0,D2y,BX,BY,'spline');
    
    %% refined phase
    if strcmp(option{3}, 'numerical_phase')
        ttao = interp2(BX0,BY0,tao,BX,BY,'spline');
        mid = round(size(tao,1)/2);
        taox = tao2x ./ (2*tao);   taox(mid, mid) = 0;
        taoy = tao2y ./ (2*tao);   taoy(mid, mid) = 0;
        ttaox = interp2(BX0,BY0,taox,BX,BY,'spline'); % refined phase
        ttaoy = interp2(BX0,BY0,taoy,BX,BY,'spline'); % refined phase
    end
    
end


%% Assembling
for p = 1:nQuad
    % quadrature points in the x-y coordinate
    pxy = lambda(p,1)*node(elem(:,1),:) ...
        + lambda(p,2)*node(elem(:,2),:) ...
        + lambda(p,3)*node(elem(:,3),:);
    reppxy = repmat(pxy,Nray,1);
    sxy = s_xy(reppxy(:,1),reppxy(:,2));
    
    
    % homogeneous case
    if ~iscell(option) && strcmp(option, 'homogeneous')
        x = (pxy(:,1)-xs);  y = (pxy(:,2)-ys);
        r = sqrt(x.^2 + y.^2);
        
        ub = 1i/4*besselh(0,1,omega*r);                  % Babich expression
        ub_g1 = -1i/4*omega*besselh(1,1,omega*r)./r.*x;  % partial derivative wrt x
        ub_g2 = -1i/4*omega*besselh(1,1,omega*r)./r.*y;  % partial derivative wrt y
    end
    
    % gravity case
    if iscell(option) && strcmp(option{1}, 'gravity')
        alpha = option{2};   E = option{3};
        trg = pxy';  src = [xs;ys];
        ub = lhelmfs(trg,src,alpha,E);
        
        trg = repmat([xs;ys], 1, size(pxy,1)); src = pxy';
        [~, ub_g1, ub_g2] = lhelmfs(trg,src,alpha,E,1);
        ub = ub(:);  ub_g1 = ub_g1(:);  ub_g2 = ub_g2(:);
    end
    
    % Babich case
    if  iscell(option) && strcmp(option{1}, 'Babich')
        
        Bu = [DD1(:), DD1x(:), DD1y(:), DD2(:), DD2x(:), DD2y(:)];
        Buint = interpolation2(Bx, By, Bu, pxy);
        amp1 = Buint(:,1);  amp1x = Buint(:,2);  amp1y = Buint(:,3);
        amp2 = Buint(:,4);  amp2x = Buint(:,5);  amp2y = Buint(:,6);
        
        if strcmp(option{3}, 'numerical_phase')
            Bu = [ttao(:), ttaox(:), ttaoy(:)];
            Buint = interpolation2(Bx, By, Bu, pxy);
            pha = Buint(:,1);  phax = Buint(:,2);  phay = Buint(:,3);
        end
        
        if strcmp(option{3}, 'exact_phase')
            if strcmp(option{2}, 'CGV')
                [~, pha, phax, phay] = eikonal_cgv(1, [0.1, -0.2], [0,0], pxy);
            end
            
            if strcmp(option{2}, 'Homo')
                x = (pxy(:,1)-xs);  y = (pxy(:,2)-ys);
                pha = sqrt(x.^2 + y.^2);
                phax = x./pha;  phay = y./pha;
            end
        end
                    
        
        c1 = 1i*(sqrt(pi)/2);
        f1 = c1*besselh(0,1,omega*pha);
        f1x = - c1*besselh(1,1,omega*pha)*omega.*phax;
        f1y = - c1*besselh(1,1,omega*pha)*omega.*phay;
        
        G1 = f1.*amp1;
        G1x = f1x.*amp1 + f1.*amp1x;
        G1y = f1y.*amp1 + f1.*amp1y;
        
        c2 = 1i*(sqrt(pi)/2)*exp(1i*pi);
        f2 = c2*(2*pha/omega).*besselh(1,1,omega*pha);
        temp = c2* ( 4/omega*besselh(1,1,omega*pha) - 2*pha.*besselh(2,1,omega*pha) );
        f2x = temp.*phax;
        f2y = temp.*phay;
        
        G2 = f2.*amp2;
        G2x = f2x.*amp2 + f2.*amp2x;
        G2y = f2y.*amp2 + f2.*amp2y;
        
        ub = G1 + G2;
        ub_g1 = G1x + G2x;
        ub_g2 = G1y + G2y;
        
    end
    
    
    
    fpxy = singularity_RHS(epsilon,xs,ys,pxy,ub,ub_g1,ub_g2);
    fp = repmat(fpxy,Nray,1).*sxy;
    
    
    for i = 1:3
        gradtempi = - ray(elem(:,i),:);
        gradtempi = gradtempi(:);
        gradtempi = [real(gradtempi), imag(gradtempi)];
        fphasei = gradtempi(:,1).*reppxy(:,1) + gradtempi(:,2).*reppxy(:,2);
        kki = kk(elem(:,i),:);
        kki = kki(:);
        phasei = exp(1i*kki.*fphasei);
        bt(:,i) = bt(:,i) + weight(p)*phi(p,i)*phasei.*fp;
    end
end

bt = bt.*repmat(reparea,1,3);
for ni = 1:Nray
    ii = (ni-1)*NT+1:1:ni*NT;
    jj = (ni-1)*N+1:1:ni*N;
    btii = bt(ii,:);
    b(jj) = accumarray(elem(:),btii(:),[N 1]);
end
clear area bt btii elem fp fphasei gradtempi ii jj k kk kki;
clear node phasei pcy pxy ray reparea reppxy;


