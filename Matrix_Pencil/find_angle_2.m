function angle = find_angle_2(angle1, angle2, tol)
%% angle1 and angle2 are two 2*Nray x 1 vector

angles = [angle1(:);angle2(:)];
angles = sort(angles);
diff_ang = angles(2:end) - angles(1:end-1);
index = find(diff_ang < tol);
angle = (angles(index+1) + angles(index))/2;
angle = angle(:); 
