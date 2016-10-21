% tic;
% A1 = assemble_Helmholtz_matrix_with_multiple_rays_original(node,elem,omega,wpml,sigmaMax,speed,ray,fquadorder);
% toc;
% tic;
% A2 = assemble_Helmholtz_matrix_with_multiple_rays_optimized(node,elem,omega,wpml,sigmaMax,speed,ray,fquadorder);
% toc;



ch = 1/100;
[cnode,celem] = squaremesh([-1,1,-1,1],ch);
cV = sin(cnode(:,1)).*cos(cnode(:,1));

h = 1/1000;
[node,elem] = squaremesh([-1,1,-1,1],h);

tic;
int1 = interpolation(cnode,celem,node,cV);
toc;

[X,Y] = meshgrid(-1:ch:1);
V = sin(X).*cos(Y);

[Xq,Yq] = meshgrid(-1:h:1);
tic;
int2 = interp2(X,Y,V,Xq,Yq);
toc;