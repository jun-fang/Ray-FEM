function fnodes = nodes_in_element(nt,node,elem,k)
fnodes = [];
if k == 1
    NT = size(elem,1);
    [node,elem] = uniformrefine(node,elem);
    fnodes = elem(nt+3*NT,:);
    fnodes = fnodes(:);
end

if k>1
    NT = size(elem,1);
    nt = [nt;nt+NT;nt+2*NT;nt+3*NT];
    [node,elem] = uniformrefine(node,elem);
    fnodes = [fnodes;nodes_in_element(nt,node,elem,k-1)];
end

