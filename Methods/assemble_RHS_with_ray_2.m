function b = assemble_RHS_with_ray_2(node,elem,omega,source,speed,ray,fquadorder)
%% Function to assemble the right hand side : 
%         -\Delta u - (omega/c)^2 u = f               in D
%                                 u = 0               on \partial D 
%  Version 2: The numbers of rays crossing at each grid node could be different!!! 
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


%% FEM set up
N = size(node,1);         % number of grid nodes
NT = size(elem,1);        % number of triangle elements
c = speed(node);      % medium speed 
k = omega./c;             % wavenumber

%% Numerical Quadrature
[lambda,weight] = quadpts(fquadorder);
phi = lambda;             % the value of piesewise linear basis functions at numerical quadrature points
nQuad = size(lambda,1);


%% Compute geometric quantities and gradient of local basis
[~,area] = gradbasis(node,elem);


%% ray information
ray_num = zeros(N,1);     % number of rays at each grid point
ray_dof = zeros(N,1);     % ray_dof(n) = sum(ray_num(1:n)) 

temp = 0;
for n = 1:N
    ray_num(n) = size(ray{n},2);
    ray_dof(n) = temp + ray_num(n);
    temp = ray_dof(n);
end

Ndof = temp;              % total degree of freedom


%% Assembling the right hand side b

b = zeros(Ndof,1);   

for nt = 1:NT    
    for p = 1:nQuad
        % quadrature points in the x-y coordinate
        pxy = lambda(p,1)*node(elem(nt,1),:) ...
            + lambda(p,2)*node(elem(nt,2),:) ...
            + lambda(p,3)*node(elem(nt,3),:);
        fp = source(pxy);
        
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

