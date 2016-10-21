function angle = find_angle_3(angle1, angle2, angle3,tol)
%% tol: angle error tolerance
ang = cell(3,1);
ang{1} = find_angle_2(angle1, angle2,tol);
ang{2} = find_angle_2(angle2, angle3,tol);
ang{3} = find_angle_2(angle3, angle1,tol);

n12 = length(ang{1});
n23 = length(ang{2});
n31 = length(ang{3});

[~,ii] = sort([n12,n23,n31]);
angle = ang{ii(end)};


