function [fb]=BGFiltrage(fu,kr,imp,L,gau,M)
%% NMLA filtering in the fourier space
% fu: FFT of impedance quantity U
% kr: k*r
% imp: parameter in impedance quantity
% L: truncation level
% gau: parameter in gaussian kernel
% M: number of sample points on the observation cicle

%% bessel and derivative of bessel
LP = max( L+2, 3) ;
Bj = besselj(0:LP-1,kr) ;                     %% bessel J_l(kr)
DBj(1) = -Bj(2) ;
DBj(2:LP-1) = 0.5*(Bj(1:LP-2) - Bj(3:LP)) ;   %% derivative of bessel

%% gausian kernel
G(1) = 1 ; 
A = gau/L ;
ii = 2:L+1 ;
G(ii) =  exp(-0.5*(A*(ii-1)).^2) ;
G = G/(2*sum(G)-1) ;

%% filtering operator
Fltr(1)  = (Bj(1)-1i*DBj(1)*imp) ;
Fltr(ii) = (Bj(ii)-1i*DBj(ii)*imp).*(1i.^(ii-1)) ;  %% ok everything shifted by 1
Fltr = G./Fltr ;

fb=zeros(1,M)+1i*0  ;
fb(1) = Fltr(1)*fu(1) ;                       %% FU_0
fb(ii) = Fltr(ii).*fu(ii)  ;                  %% FU_{1,...,L} 
fb(M-ii+2) = Fltr(ii).*fu(M-ii+2) ;           %% FU_{-1,...,-L} 

