function [node,elem] = squaremesh_annulus(square1,square2,h)
%% SQUAREMESH uniform mesh of a square
%
% [node,elem] = squaremesh([x0,x1,y0,y1],h) generates a uniform mesh of the
% rectangle [x0,x1]*[y0,y1] with mesh size h.
%
% Example
%
%   [node,elem] = squaremesh([0,1,0,1],0.2);
%   showmesh(node,elem);
%   findnode(node);
%   findelem(node,elem);
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

%% Generate nodes
x0 = square1(1); x1 = square1(2); 
y0 = square1(3); y1 = square1(4);
[x,y] = meshgrid(x0:h:x1,y0:h:y1);
idx = (x> square2(1)+10*eps).*(x< square2(2) - 10*eps).*...
    (y> square2(1)+10*eps).*(y< square2(4) - 10*eps);
idx = 1-idx;
idx = idx(:);
idx = find(idx);
x1 = x(:); y1 = y(:);
node1 = [x1,y1];
node = [x1(idx), y1(idx)];


%% Generate elements
ni = size(x,1); % number of rows
N = length(x1);
t2nidxMap = 1:N-ni;
topNode = ni:ni:N-ni;
t2nidxMap(topNode) = [];
k = (t2nidxMap)';
nodek = node1(k,:);
xk = nodek(:,1); yk = nodek(:,2);
idxk = (xk>= square2(1)-10*eps).*(xk< square2(2) - 10*eps).*...
    (yk>= square2(1)- 10*eps).*(yk< square2(4) - 10*eps);
idxk = 1-idxk;
k = k(idxk>0.5);
kl = k + ni;

[~,k,~] = intersect(idx,k);
[~,kl,~] = intersect(idx,kl);

elem = [kl kl+1 k; k+1 k kl+1];
% 4 k+1 --- k+ni+1 3  
%    |        |
% 1  k  ---  k+ni  2