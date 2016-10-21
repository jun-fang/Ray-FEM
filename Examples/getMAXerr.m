function err = getMAXerr(u,h,ru,rh,wpml)
%% Maximum error of the approximation u to the reference solution uh on square physical domain
N = length(u);
n = round(sqrt(N));
u = reshape(u,n,n);

rN = length(ru);
rn = round(sqrt(rN));
ru = reshape(ru,rn,rn);

r = round(h/rh);

idx = 1:n;
idx = (idx-1)*r+1;

du = u - ru(idx,idx);
nwpml = ceil(wpml/h);
phy_du = du(nwpml+1:n-nwpml, nwpml+1:n-nwpml);
err = norm(phy_du(:),inf);