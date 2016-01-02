function diff = angle_error(angle1,angle2)
d1 = angle1 - angle2;
d2 = d1 + 2*pi;
d3 = d1 - 2*pi;
diff = min(abs(d1),abs(d2));
diff = min(diff,abs(d3));