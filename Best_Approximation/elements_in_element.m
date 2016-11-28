function [rnt] = elements_in_element(nt,node,elem,k)

NT = size(elem,1);
rnt = [nt;nt+NT;nt+2*NT;nt+3*NT];

if k>1
    NT = size(elem,1);
    [node,elem] = uniformrefine(node,elem);
    nt = [nt;nt+NT;nt+2*NT;nt+3*NT];
    rnt = elements_in_element(nt,node,elem,k-1);
end    
    