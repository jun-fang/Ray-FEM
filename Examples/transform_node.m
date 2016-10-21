function [new_node,new_u] = transform_node(node,u)
N = length(node(:,1));
n = round(sqrt(N));
[X,I] = sort(node(:,1));
Y = node(I,2);
u = u(I);
for i = 1:n
    [y,ii] = sort(Y(((i-1)*n+1):(i*n)));
    Y(((i-1)*n+1):(i*n)) = y;
    u(((i-1)*n+1):(i*n)) = u((i-1)*n + ii);
end
new_node = [X,Y];
new_u = u;
    