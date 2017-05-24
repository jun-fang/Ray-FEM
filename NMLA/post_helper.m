function idx = post_helper(angles)
%% helper function to indentify which qudrant the angles locate

% angles: 1xNray, Nray <= 4
% 0 <= angles <= 2*pi

b1 = angles <= pi/2;
b2 = (pi/2 <= angles) .* (angles < pi);
b3 = (pi <= angles) .* (angles < 3*pi/2);
b4 = 1 - (b1 + b2 + b3);
idx = 1*b1 + 2*b2 + 3*b3 + 4*b4;
n = length(idx);
if n > 1
    for i = 1:n-1
        if idx(i) == idx(i+1) 
            if idx(i) < 4
                idx(i+1) = idx(i+1) + 1;
            else
                idx(i) = idx(i) - (n-i);
                idx(i+1) = idx(i+1) - (n-i-1);
            end
        end
    end
end

if idx(end) > 4
    idx = idx - (idx-4);
end