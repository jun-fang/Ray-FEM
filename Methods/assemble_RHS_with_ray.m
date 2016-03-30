function b = assemble_RHS_with_ray(node,elem,omega,source,speed,ray,fquadorder)
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
fprintf('Assembling the right-hand side \n');

%% FEM set up
N = size(node,1);        % number of grid points
NT = size(elem,1);       % number of triangle elements
Nray = size(ray,2);
Ndof = N*Nray;          % degree of freedom 

k = omega./speed(node);          % wavenumber
kk = repmat(k,1,Nray);

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
    
    % we suppose that the source is well inside the physical domain
    fp = source(pxy);
%     fp = source(pxy).*( pxy(:,1) < xmax - wpml ).*( pxy(:,1) > xmin + wpml )...
%         .*( pxy(:,2) < ymax - wpml ).*( pxy(:,2) > ymin + wpml ); 
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
clear area b bt btii elem fp fphasei gradtempi ii jj k kk kki;
clear node phasei pcy pxy ray reparea reppxy;

frpintf('\n');

