function [node,elem] = symmetric_squaremesh(square,n)
%% symmetric mesh of a square 

% [node,elem] = symmetric_squaremesh([x0,x1,y0,y1],n) generates a mesh of the
% rectangle [x0,x1]*[y0,y1] with mesh size hx = (x1-x0)/2^n, hy = (y1-y0)/2^n.

x0 = square(1); x1 = square(2); 
y0 = square(3); y1 = square(4);
node = [x0,y0; x1,y0; x1,y1; x0,y1];
elem = [2,3,1; 4,1,3];
for k = 1:n
    [node,elem] = uniformbisect(node,elem);
end