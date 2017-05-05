function [ux,uy] = num_derivative(u,h,order,m,n)
%% numerical central difference schemes to compute the partial derivatives of u
% u is a N x 1 vector

nr = size(u,1);  nc = size(u,2);
N = length(u(:));
if nargin < 4
    n = round(sqrt(N)); m = n;
end

u = reshape(u,n,m);
ux = u; uy = u;

if order == 2   % convergbence order O(h^2)
    ux(:,2:m-1) = (u(:,3:m) - u(:,1:m-2))/(2*h);
    ux(:,m) = 2*ux(:,m-1) - ux(:,m-2);
    ux(:,1) = 2*ux(:,2) - ux(:,3);    
    
    uy(2:n-1,:) = (u(3:n,:) - u(1:n-2,:))/(2*h);
    uy(n,:) = 2*uy(n-1,:) - uy(n-2,:);
    uy(1,:) = 2*uy(2,:) - uy(3,:); 
end

if order == 4   % n >= 8, convergbence order O(h^4)
    ux(:,3:m-2) = ( -u(:,5:m) + 8*u(:,4:m-1) - 8*u(:,2:m-3) + u(:,1:m-4) )/(12*h);
    ux(:,m-1) = 4*ux(:,m-3) - ux(:,m-2) - ux(:,m-4) - ux(:,m-5);
    ux(:,m) = 4*ux(:,m-2) - ux(:,m-1) - ux(:,m-3) - ux(:,m-4);
    ux(:,2) = 4*ux(:,4) - ux(:,3) - ux(:,5) - ux(:,6);
    ux(:,1) = 4*ux(:,3) - ux(:,2) - ux(:,4) - ux(:,5);
    
    uy(3:n-2,:) = ( -u(5:n,:) + 8*u(4:n-1,:) - 8*u(2:n-3,:) + u(1:n-4,:) )/(12*h);
    uy(n-1,:) = 4*uy(n-3,:) - uy(n-2,:) - uy(n-4,:) - uy(n-5,:);
    uy(n,:) = 4*uy(n-2,:) - uy(n-1,:) - uy(n-3,:) - uy(n-4,:);
    uy(2,:) = 4*uy(4,:) - uy(3,:) - uy(5,:) - uy(6,:);
    uy(1,:) = 4*uy(3,:) - uy(2,:) - uy(4,:) - uy(5,:);
end

ux = reshape(ux, nr, nc);  
uy = reshape(uy, nr, nc);