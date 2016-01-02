function A = tri_area(node1,node2,node3)
%% Computing the area of triangles

vec2 = node1 - node3;
vec3 = node2 - node1;

area = 0.5*(-vec3(:,1).*vec2(:,2) + vec3(:,2).*vec2(:,1));
A = abs(area);

end
