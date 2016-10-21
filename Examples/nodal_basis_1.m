function f = nodal_basis_1(xc,yc,h,xy)
%%  Nodal basis function of the following mesh:
% 
%  xc = 0; yc = 0; h = 1;  
%  node = [xc-h,yc-h;xc+h,yc-h;xc+h,yc+h;xc-h,yc+h];
%  elem = [2,3,1; 4,1,3];
%  [node,elem] = uniformbisect(node,elem);
%  showmesh(node,elem)
%
% Example: 
%  [node,elem] = squaremesh([-1,1,-1,1],1/32);
%  f = nodal_basis_1(0,0,1/2,node);
%  FJ_showresult(node,elem,f);
% 
x = xy(:,1); y = xy(:,2);
mxy = max(abs(x-xc), abs(y-yc));
idx = mxy<h;
f = (1 - mxy/h).*idx;