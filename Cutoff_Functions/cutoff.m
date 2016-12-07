function cf = cutoff(node,xs,ys,a,b)
%% Smooth cut-off function 
% Input: 
%     node: Nx2 matrix, x=node(:,1) and y=node(:,2) are x,y coordinates
%     (xs,ys): center, r = sqrt((node(:,1)-xs).^2 + (node(:,2)-yx).^2)
%     a,b: positve numbers such that 0<a<b
% 
% Ouput:  
%     cf: Nx1 vector, 
%           
%     cf(j) = 1        if r<=a
%             smooth   if a<r<b
%             0        if r>=b

x = node(:,1);  y = node(:,2);
r = sqrt((x-xs).^2 + (y-ys).^2);
t = (r-a)/(b-a);
cf = basic_cutoff(t);
cf(r<=a) = 1;   cf(r>=b) = 0;
