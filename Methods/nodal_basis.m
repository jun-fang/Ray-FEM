function f = nodal_basis(xc,yc,h,xy)
%%  Nodal basis function of the following mesh:
% 
%  xc = 0; yc = 0; h = 1;  
%  node = [xc,yc; xc-h,yc; xc,yc+h; xc+h,yc+h; xc+h,yc; xc,yc-h; xc-h,yc-h;];
%  elem = [1,2,3; 3,4,1; 1,4,5; 5,6,1; 1,6,7; 7,2,1];
%  showmesh(node,elem)
%
% Example: 
%  [node,elem] = squaremesh([-1,1,-1,1],1/32);
%  f = nodal_basis(0,0,1/2,node);
%  FJ_showresult(node,elem,f);
% 

x = xy(:,1); y = xy(:,2);
xx = x - xc; yy = y - yc;
ax = abs(xx); ay = abs(yy);
idx = xx.*yy>0;
f = (1 - max(ax,ay)/h).*idx.*(max(ax,ay)<h) ...
    + (1-(ax+ay)/h).*(1-idx).*((ax+ay)<h);
