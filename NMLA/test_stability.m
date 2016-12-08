
sta_coef = 0.775;

nt = 1000;
krs = zeros(1,nt);
errs = zeros(1,nt);
for n = 1:nt
    kr = 10*n;
    krs(n) = kr;
    L = floor(sta_coef*sqrt(kr));
    idx = -L:L;
    Cl = (1i).^idx.*besselh(idx,1,kr);
    Cl = sqrt(1i*pi*kr/2)*exp(-1i*kr)*Cl;
    err = sum(abs(Cl-1))/(2*L+1);
    errs(n) = err;
end

plot(krs,errs)