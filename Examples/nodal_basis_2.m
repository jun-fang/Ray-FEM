function f = nodal_basis_2(xc,yc,h,xy)
%%  Nodal basis function of the following mesh:
% 
%  xc = 0; yc = 0; h = 1;  
%  node = [xc,yc; xc-h,yc; xc,yc+h; xc+h,yc; xc,yc-h];
%  elem = [1,2,3; 3,4,1; 1,4,5; 5,2,1];
%  showmesh(node,elem)
%
% Example: 
%  [node,elem] = squaremesh([-1,1,-1,1],1/32);
%  f = nodal_basis_2(0,0,1/2,node);
%  FJ_showresult(node,elem,f);
%
x = xy(:,1); y = xy(:,2);
xx = abs(x-xc);
yy = abs(y-yc);
xxyy = xx+yy;
idx = xxyy<h;
f = (1 - xxyy/h).*idx;