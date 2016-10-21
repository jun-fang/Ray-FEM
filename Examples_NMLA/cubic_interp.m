omega = 20*pi:20*pi:400*pi;
omega = 10*pi:10*pi:80*pi;
n = length(omega);
Rest = 1;
r = zeros(1,n);
for i = 1:n
%     p = [1,0,1,-1.5-0.56*(omega(i)*Rest)^0.75];
p = [1,0,1,-2.5-0.93*(omega(i)*Rest)^0.763];
rt = roots(p);
pidx = find(rt>0);
r(i) = (rt(pidx(1)))^3/omega(i);  
end

showrate(omega,1./(omega.*r))


for i = 1:n
    Rest = 100 + omega(i)/2;
%     p = [1,0,1,-1.5-0.56*(omega(i)*Rest)^0.75];
p = [1,0,1,-2.5-0.93*(omega(i)*Rest)^0.763];
rt = roots(p);
pidx = find(rt>0);
r(i) = (rt(pidx(1)))^3/omega(i);  
end

showrate(omega,1./(omega.*r))
