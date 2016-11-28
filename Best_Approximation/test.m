[node,elem] = squaremesh([0,1,0,1],1);
subplot(2,2,1);
findelem(node,elem);
findnode(node);
nodes_in_element(1,node,elem,2)
[node,elem] = uniformrefine(node,elem);
subplot(2,2,2);

findelem(node,elem);
findnode(node);

[node,elem] = uniformrefine(node,elem);
% [node,elem] = squaremesh([0,1,0,1],1/2);
subplot(2,2,3);
findelem(node,elem);
findnode(node);

[node,elem] = uniformrefine(node,elem);
subplot(2,2,4);

findelem(node,elem);
findnode(node);
