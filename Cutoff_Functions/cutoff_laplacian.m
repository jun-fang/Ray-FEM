function l = cutoff_laplacian(a,b,node,xs,ys)
%% Laplacian of smooth cut-off function 
% Input: 
%     a,b: positve numbers such that 0<a<b
%     node: Nx2 matrix, x=node(:,1) and y=node(:,2) are x,y coordinates
%     (xs,ys): center, r = sqrt((node(:,1)-xs).^2 + (node(:,2)-yx).^2)
% 
% Output:  
%     l: Nx1 vector, 
%           
%     l(j) = [0,0]      if r<=a or r>=b
%              smooth     if a<r<b

x = node(:,1);   y = node(:,2);
r = sqrt((x-xs).^2 + (y-ys).^2);
t = (r-a)/(b-a);
s1 = basic_cutoff_first_derivative(t);
s2 = basic_cutoff_second_derivative(t);

l = s2/(b-a)^2 + s1./((b-a)*r);
l(r<=a) = 0;   l(r>=b) = 0;