function err = getRelL2err(u,h,ru,rh,wpml)
%% Relative L2 error of the approximation u to the reference solution uh on square physical domain
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

% cc = (n-1)/2 + 1;
% du(cc-16:cc+16,cc-16:cc+16) = 0;

nwpml = ceil(wpml/h);
rnwpml = ceil(wpml/rh);
phy_du = du(nwpml+1:n-nwpml, nwpml+1:n-nwpml);
phy_ru = ru(rnwpml+1:rn-rnwpml, rnwpml+1:rn-rnwpml);
err = (norm(phy_du(:),2)*h)/(norm(phy_ru(:),2)*rh);