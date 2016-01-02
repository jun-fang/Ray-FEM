function [ux,uy] = num_derivative_2(u,h,m,n,order)
%% Numerical central difference schemes to compute the derivative of u
% u: N x 1 vector, N = m*n, m and n are the numbers of grid points on each
%    dimension
% h: mesh size
% order: convergence order of numerical scheme
%
% ux,uy: N x 1 vectors, partial derivatives
%        \partial u / \partial x, \partial u / \partial y 

u = reshape(u,m,n);
ux = u; uy = u;

if order == 2
    uy(2:m-1,:) = (u(3:m,:) - u(1:m-2,:))/(2*h);
    uy(m,:) = 2*uy(m-1,:) - uy(m-2,:);
    uy(1,:) = 2*uy(2,:) - uy(3,:); 
    
    ux(:,2:n-1) = (u(:,3:n) - u(:,1:n-2))/(2*h);
    ux(:,n) = 2*ux(:,n-1) - ux(:,n-2);
    ux(:,1) = 2*ux(:,2) - ux(:,3);        
end

if order == 4   % n >= 8
    uy(3:m-2,:) = ( -u(5:m,:) + 8*u(4:m-1,:) - 8*u(2:m-3,:) + u(1:m-4,:) )/(12*h);
    uy(m-1,:) = 4*uy(m-3,:) - uy(m-2,:) - uy(m-4,:) - uy(m-5,:);
    uy(m,:) = 4*uy(m-2,:) - uy(m-1,:) - uy(m-3,:) - uy(m-4,:);
    uy(2,:) = 4*uy(4,:) - uy(3,:) - uy(5,:) - uy(6,:);
    uy(1,:) = 4*uy(3,:) - uy(2,:) - uy(4,:) - uy(5,:);
    
    ux(:,3:n-2) = ( -u(:,5:n) + 8*u(:,4:n-1) - 8*u(:,2:n-3) + u(:,1:n-4) )/(12*h);
    ux(:,n-1) = 4*ux(:,n-3) - ux(:,n-2) - ux(:,n-4) - ux(:,n-5);
    ux(:,n) = 4*ux(:,n-2) - ux(:,n-1) - ux(:,n-3) - ux(:,n-4);
    ux(:,2) = 4*ux(:,4) - ux(:,3) - ux(:,5) - ux(:,6);
    ux(:,1) = 4*ux(:,3) - ux(:,2) - ux(:,4) - ux(:,5);
    
end

ux = ux(:);  uy = uy(:);