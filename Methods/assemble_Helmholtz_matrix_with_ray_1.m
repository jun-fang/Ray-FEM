function [A] = assemble_Helmholtz_matrix_with_ray_1(node,elem,omega,wpml,sigmaMax,speed,ray,fquadorder)
%% Function to assemble the Ray-FEM Helmholtz matrix with PML: 
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
%   ray: N x Nray matrix that contains the ray information
%        stored as the complex form exp(i*ray_angle)
% 
%   fquadorder: The order of numerical quadrature
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% OUTPUT:
%   
%   A: Ndof x Ndof Helmholtx matrix, Ndof = N*Nray
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% FEM set up
N = size(node,1);        % number of grid points
%NT = size(elem,1);       % number of triangle elements
Nray = size(ray,2);     % number of rays crossing at each grid node
Ndof = N*Nray;          % degree of freedom 

k = omega./speed(node);           % wavenumber

% computing the limits of the domain
xmax = max(node(:,1));
xmin = min(node(:,1));
ymax = max(node(:,2));
ymin = min(node(:,2));


% %% PML set up
% sigmaPML_x = @(x)sigmaMax*( (x-xmin-wpml).^2.*(x < xmin + wpml) + ...
%                 (x-(xmax-wpml)).^2.*(x > xmax - wpml))/wpml^2;         
% sigmaPML_y = @(y) sigmaMax*( (y-ymin-wpml).^2.*(y < ymin + wpml) ...
%                 + (y-(ymax-wpml)).^2.*(y > ymax - wpml))/wpml^2;
% s_x = @(x,y) (1+1i*sigmaPML_y(y)/omega)./(1+1i*sigmaPML_x(x)/omega);       %% s1/s2
% s_y = @(x,y) (1+1i*sigmaPML_x(x)/omega)./(1+1i*sigmaPML_y(y)/omega);       %% s2/s1
% s_xy = @(x,y) ((1+1i*sigmaPML_x(x)/omega).*(1+1i*sigmaPML_y(y)/omega));    %% 1/(s1*s2)

% in order to define the correct boundary conditionts at the interfaces
% we need to specify the width of the PML in eevry direction
%    ---wpml(1)---
%   |             |
%   |             |
%  wpml(3)     wpml(4)
%   |             |
%   |             |
%    ---wpml(2)---
%

% if the different lengths of the pml are not specified we specify all of
% them uisng the same number of pml points

if length(wpml) == 1
    wpml = ones(4,1)*wpml;
end

%% PML set up

% usual quadratic profile 
sigmaPML_x = @(x) sigmaMax*( ( (x-xmin-wpml(1)) / wpml(1) ).^2.*(x < xmin + wpml(1)) + ...
                ( (x-(xmax-wpml(2))) / wpml(2) ).^2.*(x > xmax - wpml(2)) );         
sigmaPML_y = @(y) sigmaMax*( ( (y-ymin-wpml(3)) / wpml(3) ).^2.*(y < ymin + wpml(3)) + ...
                ( (y-(ymax-wpml(4))) / wpml(4) ).^2.*(y > ymax - wpml(4) ) );

s_x = @(p) (1+1i*sigmaPML_y(p(:,2))/omega)./(1+1i*sigmaPML_x(p(:,1))/omega);       %% s1/s2
s_y = @(p) (1+1i*sigmaPML_x(p(:,1))/omega)./(1+1i*sigmaPML_y(p(:,2))/omega);       %% s2/s1
s_xy = @(p) ((1+1i*sigmaPML_x(p(:,1))/omega).*(1+1i*sigmaPML_y(p(:,2))/omega));    %% 1/(s1*s2)

% % unbounded PML profile             
% sigmaPML_x = @(p) speed(p).*( 1./( p(:,1)-xmin + 10*eps ).*(p(:,1) < xmin + wpml(1) + 10*eps) + ...
%                 1./( xmax-p(:,1) + 10*eps ).*( p(:,1) > xmax - wpml(2) - 10*eps) ) ;         
% sigmaPML_y = @(p) speed(p).*( 1./( p(:,2)-ymin + 10*eps ).*(p(:,2) < ymin + wpml(3) + 10*eps) + ...
%                 1./( ymax-p(:,2) + 10*eps ).*( p(:,2) > ymax - wpml(4) - 10*eps) );
%               
%             
% s_x = @(p) (1+1i*sigmaPML_y(p)/omega)./(1+1i*sigmaPML_x(p)/omega);       %% s1/s2
% s_y = @(p) (1+1i*sigmaPML_x(p)/omega)./(1+1i*sigmaPML_y(p)/omega);       %% s2/s1
% s_xy = @(p) ((1+1i*sigmaPML_x(p)/omega).*(1+1i*sigmaPML_y(p)/omega));    %% 1/(s1*s2)


%% Numerical Quadrature
[lambda,weight] = quadpts(fquadorder);
phi = lambda;           % linear bases
nQuad = size(lambda,1);


%% Compute geometric quantities and gradient of local basis
[Dphi,area] = gradbasis(node,elem);


%% Assembling the matrix

A = sparse(Ndof,Ndof);
% M = sparse(Ndof,Ndof);

for p = 1:nQuad
    % quadrature points in the x-y coordinate
    pxy = lambda(p,1)*node(elem(:,1),:) ...
        + lambda(p,2)*node(elem(:,2),:) ...
        + lambda(p,3)*node(elem(:,3),:);
    
    % building the PML profiles
    sx = s_x(pxy);
    sy = s_y(pxy);
    sxy = s_xy(pxy);
    
    % local wavenumber
    k2 = (omega./speed(pxy)).^2;
    
    for i = 1:3   
        for j = 1:3
            for nii = 1: Nray
                % phase e^{-ik ray_dierection \dot pxy}
                gradtempi = - ray(elem(:,i),nii);
                gradtempi = [real(gradtempi), imag(gradtempi)];
                fphasei = gradtempi(:,1).*pxy(:,1) + gradtempi(:,2).*pxy(:,2);
                ki = k(elem(:,i));
                phasei = exp(1i*ki.*fphasei);
                
                for njj = 1: Nray
                    % phase e^{ik ray_dierection \dot pxy}
                    gradtempj = ray(elem(:,j),njj);
                    gradtempj = [real(gradtempj), imag(gradtempj)];
                    fphasej = gradtempj(:,1).*pxy(:,1) + gradtempj(:,2).*pxy(:,2);
                    kj = k(elem(:,j));
                    phasej = exp(1i*kj.*fphasej);
                    
                    exp_phase = phasei.*phasej;
                    
                    tempA1 = sx.*Dphi(:,1,i).*Dphi(:,1,j) + sy.*Dphi(:,2,i).*Dphi(:,2,j); 
                    
                    tempA2 = 1i*ki*phi(p,i).*(sx.*gradtempi(:,1).*Dphi(:,1,j) + sy.*gradtempi(:,2).*Dphi(:,2,j))...
                        + 1i*kj*phi(p,j).*(sx.*gradtempj(:,1).*Dphi(:,1,i) + sy.*gradtempj(:,2).*Dphi(:,2,i));  
                    
                    tempA3 = phi(p,i)*phi(p,j)*ki.*kj.*(sx.*gradtempi(:,1).*gradtempj(:,1)...
                        + sy.*gradtempi(:,2).*gradtempj(:,2));
                    
                    tempM = sxy*phi(p,i)*phi(p,j);
                    
                    rows = (nii-1)*N + elem(:,i);
                    cols = (njj-1)*N + elem(:,j);
                    vals = weight(p)*(tempA1 + tempA2 - tempA3 - k2.*tempM).*exp_phase.*area;
                    
                    A = A + sparse(rows,cols,vals,Ndof,Ndof);
%                     M = M + sparse(rows,cols,weight(p)*(tempM).*exp_phase.*area,Ndof,Ndof);
                end
            end
        end
    end
end

clear Dphi fphasei fphasej phasei phasej exp_phase gradtempi gradtempj;
clear k2 ki kj pxy tempA1 tempA2 tempA3 tempM;
