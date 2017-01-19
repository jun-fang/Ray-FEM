function g = cutoff_gradient(a,b,node,xs,ys)
%% Gradient of smooth cut-off function 
% Input: 
%     a,b: positve numbers such that 0<a<b
%     node: Nx2 matrix, x=node(:,1) and y=node(:,2) are x,y coordinates
%     (xs,ys): center, r = sqrt((node(:,1)-xs).^2 + (node(:,2)-yx).^2)
% 
% Output:  
%     g: Nx2 matrix, 
%           
%     g(j,:) = [0,0]      if r<=a or r>=b
%              smooth     if a<r<b

x = node(:,1);   y = node(:,2);
r = sqrt((x-xs).^2 + (y-ys).^2);
rx = (x-xs)./r;  ry = (y-ys)./r;
t = (r-a)/(b-a);
s1 = basic_cutoff_first_derivative(t);
g1 = s1.*rx/(b-a);  g1(r<=a) = 0;   g1(r>=b) = 0;
g2 = s1.*ry/(b-a);  g2(r<=a) = 0;   g2(r>=b) = 0;
g = [g1, g2];
