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

for p = 1:nQuad
    % quadrature points in the x-y coordinate
    pxy = lambda(p,1)*node(elem(:,1),:) ...
        + lambda(p,2)*node(elem(:,2),:) ...
        + lambda(p,3)*node(elem(:,3),:);
    reppxy = repmat(pxy,Nray,1);
    sxy = s_xy(reppxy(:,1),reppxy(:,2));
    
    
    % homogeneous case
    if option == 'homogeneous'
        x = (pxy(:,1)-xs);  y = (pxy(:,2)-ys);
        r = sqrt(x.^2 + y.^2);
        
        ub = 1i/4*besselh(0,1,omega*r);                  % Babich expression
        ub_g1 = -1i/4*omega*besselh(1,1,omega*r)./r.*x;  % partial derivative wrt x
        ub_g2 = -1i/4*omega*besselh(1,1,omega*r)./r.*y;  % partial derivative wrt y
    end
    
    % gravity case
    if iscell(option) && option{1} == 'gravity'
        alpha = option{2};   E = option{3};
        trg = pxy';  src = [xs;ys];
        ub = lhelmfs(trg,src,alpha,E);
        
        trg = repmat([xs;ys], 1, size(pxy,1)); src = pxy';
        [~, ub_g1, ub_g2] = lhelmfs(trg,src,alpha,E,1);
        ub = ub(:);  ub_g1 = ub_g1(:);  ub_g2 = ub_g2(:);
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


