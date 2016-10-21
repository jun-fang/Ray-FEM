function r = NMLA_radius(omega,Rest)
r = 0*omega;
for i = 1:length(omega)
    poly = [1,0,1,-1.5-0.56*(omega(i)*Rest).^0.75];
%     poly = [1,0,1,-1.5-0.775*(omega(i)*Rest)^0.5];
%     poly = [1,0,1,-1.5-0.85*(omega(i)*Rest)^0.5];
    rt = roots(poly);
    pidx = find(rt>0);
    r(i) = rt(pidx(1))^3/omega(i);
end