% a = 1/2;
% bd = 1/4;
% h = 1/16;
% %[node,elem] = squaremesh([-a,a,-a,a],h);
% [node,elem] = squaremesh_annulus([-a,a,-a,a],[-bd,bd,-bd,bd],h);


x = 0;
y = 0;
r = 1;
h = 0.1;

[node,elem] = annulusmesh(0,0,1,0.5,0.05);
%[node,elem] = circlemesh(0,0,1,0.05);
showmesh(node,elem);
