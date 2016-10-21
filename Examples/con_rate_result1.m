omegas = [4 8 16 32]'*2*pi;
NPW = 16;
fquadorder = 9;                  % numerical quadrature order
wpml = 1/4;
sigmaMax = 25/wpml;
rh = 1/2048;
a = 1/2;
sigma = 2*h;
source = @(x) 1/(2*pi*sigma^2)*exp(-( (x(:,1)-xs).^2 + (x(:,2)-ys).^2 )/2/sigma^2);


idx = find(d<16*sigma);
    for i = 1:length(idx)
        ii = idx(i);
        xi = node(ii,1); yi = node(ii,2);
        fb(ii) = RHS_integral_with_ray(xi,yi,h,rh,omega,ray(ii),source,fquadorder);
    end

fsb = sb;
    for i = 1:length(idx)
        ii = idx(i);
        xi = node(ii,1); yi = node(ii,2);
        fsb(ii) = RHS_integral(xi,yi,h,rh,source,fquadorder);
    end


--------------------------------------------------------------------------------
Gaussian point source parameter sigma = 1/8 wavelength

--------------------------------------------------------------------------------
Ray-FEM:
 omega/(2*pi) = 04,   1/h = 64   NPW = 16 



rfem_max_err =

   0.034523404541124
   0.034347609729140
   0.034324182097472
   0.034256020311265


rfem_l2_err =

   0.145484217456133
   0.144206235630938
   0.143154704979887
   0.131421826124946


sfem_max_err =

   0.006514511493931
   0.009434968471347
   0.013096642674405
   0.017059938081681


sfem_l2_err =

   0.043071270191482
   0.081379760515742
   0.156331813527216
   0.290818381130603

   
