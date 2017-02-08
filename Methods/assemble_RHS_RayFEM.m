function b = assemble_RHS_RayFEM(node,elem,omega,wpml,sigmaMax,source,speed,ray,fquadorder)
%% Function to assemble the right hand side :
%         -\Delta u - (omega/c)^2 u = f               in D
%                                 u = 0               on \partial D
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%   omega: Angular frequency
%
%   wpml: Width of the PML
%
%   sigmaMax: Maximun absorbtion
%
%   ray: N x Nray matrix or Nx1 cell that contains the ray information
%        stored as the complex form exp(i*ray_angle)
%
%   source: function handle defining the source
%
%   fquadorder: The order of numerical quadrature
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% OUTPUT:
%
%   b: Ndof x 1 Galerking projection of the source
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fprintf('Assembling the right-hand side \n');


%% FEM set up
N = size(node,1);         % number of grid nodes
NT = size(elem,1);        % number of triangle elements
c = speed(node);      % medium speed
k = omega./c;             % wavenumber

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
phi = lambda;             % the value of piesewise linear basis functions at numerical quadrature points
nQuad = size(lambda,1);


%% Compute geometric quantities and gradient of local basis
[~,area] = gradbasis(node,elem);


%% Assembling the RHS: ray could be a N x Nray matrix or N x 1 cell

%% ray is N x Nray matrix: each grid point has the same number of rays
if ~iscell(ray)
    
    Nray = size(ray,2);
    Ndof = N*Nray;          % degree of freedom
    kk = repmat(k,1,Nray);
    reparea = repmat(area,Nray,1);
    
    % Assemble right-hand side
    bt = zeros(NT*Nray,3);       % the right hand side
    b = zeros(Ndof,1);
    
    for p = 1:nQuad
        % quadrature points in the x-y coordinate
        pxy = lambda(p,1)*node(elem(:,1),:) ...
            + lambda(p,2)*node(elem(:,2),:) ...
            + lambda(p,3)*node(elem(:,3),:);
        reppxy = repmat(pxy,Nray,1);
        sxy = s_xy(reppxy(:,1),reppxy(:,2));
        
        % we suppose that the source is well inside the physical domain
        fp = source(pxy).*sxy;
        
        if norm(fp)
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
    end
    
    bt = bt.*repmat(reparea,1,3);
    for ni = 1:Nray
        ii = (ni-1)*NT+1:1:ni*NT;
        jj = (ni-1)*N+1:1:ni*N;
        btii = bt(ii,:);
        b(jj) = accumarray(elem(:),btii(:),[N 1]);
    end
    
   
else
    %% ray is N x 1 cell: each grid point may have different number of rays
    
    ray_num = zeros(N,1);     % number of rays at each grid point
    ray_dof = zeros(N,1);     % ray_dof(n) = sum(ray_num(1:n))
    
    temp = 0;
    for n = 1:N
        ray_num(n) = size(ray{n},2);
        ray_dof(n) = temp + ray_num(n);
        temp = ray_dof(n);
    end
    
    Ndof = temp;              % total degree of freedom

    % Assemble right-hand side
    b = zeros(Ndof,1);
    
    for nt = 1:NT
        for p = 1:nQuad
            % quadrature points in the x-y coordinate
            pxy = lambda(p,1)*node(elem(nt,1),:) ...
                + lambda(p,2)*node(elem(nt,2),:) ...
                + lambda(p,3)*node(elem(nt,3),:);
            sxy = s_xy(pxy(:,1),pxy(:,2));
            fp = source(pxy).*sxy;
            
            if norm(fp)
                for i = 1:3
                    ei = elem(nt,i);    % the ei'th grid point corresponding to the i'th vertix of element nt
                    ni = ray_num(ei);   % number of rays crossing at grid point ei
                    for nii = 1: ni
                        gradtempi = - ray{ei}(nii);
                        gradtempi = [real(gradtempi), imag(gradtempi)];
                        fphasei = gradtempi(1)*pxy(1) + gradtempi(2)*pxy(2);
                        ki = k(ei);
                        phasei = exp(1i*ki*fphasei);
                        ii = ray_dof(ei) - ni + nii;
                        b(ii) = b(ii) + weight(p)*phi(p,i)*phasei*fp*area(nt);   % right hand side
                    end
                end
            end
        end
    end
    
    
end




