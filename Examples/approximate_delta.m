function f = approximate_delta(xc,yc,h,xy)
x = xy(:,1); y = xy(:,2);
mxy = max(abs(x-xc), abs(y-yc));
idx = mxy<h;
f = (1 - mxy/h).*idx;