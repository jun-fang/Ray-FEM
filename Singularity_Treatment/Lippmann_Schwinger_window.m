function window = Lippmann_Schwinger_window(y, alpha, beta)
% alpha > beta??

x = abs(y);
window = exp( 2*exp( -(alpha- beta)./( x-beta) )...
    ./ ( (x-beta)./(alpha- beta) -1 ) );

window(x <= beta+eps) = 1;
window(x >= alpha-eps) = 0;
