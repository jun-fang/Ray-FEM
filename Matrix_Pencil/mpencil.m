function [p,s] = mpencil(y,L,K)
% MP Matrix Pencil method
%
% [p,s] = mp(y,L,K),
%
% p : discrete poles (for the forward model)
% s : singular values taken from SVD
% y : time series
% L : actual order (default: autodetect, please see code)
% K : overspecified order (default: floor(N/3))
%
% example:
% t = 0:24 ;
% y = sin(3*t+pi/4)+randn(size(t))/10 ;
% z = mp(y,2,8) ;
% log(z) % the poles are z = exp(-\lambda)

% see Hua, Yingbo and Sarkar, Tapan K., 
% "Matrix Pencil Method for Estimating Parameters of 
% Exponentially Damped/Undamped Sinusoids in Noise",
% IEEE Trans. Acoust., Speech, Signal Processing, 
% Vol. ASSP-38, pp. 814-824, May 1990.

% Copyright (C) 1994 Thomas Philip Runarsson (e-mail: tpr@verk.hi.is)
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software Foundation,
% Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

% check arguments
  N = length(y) ;
  if prod(size(y)) ~= N, error('MP: y must be a single time series') ; end
  y = y(:) ;
  
% default overspecified order
  if nargin < 3, K = floor(N/3) ; end
  
% bLP data Matrix (backward linear predictor)
  for i=1:K+1, A(:,i) = y(i:i+N-K-1) ; end

% perform SVD
  [u,s,v] = svd(A(:,2:K+1),0) ;
  
% estimate true model order (you may want to plot sdd) 
% this method is not very good for high signal to noise ratios
  if nargin < 2,
    sdiag = diag(s) ; sdiag = sdiag/sdiag(1) ;
    sdd = -diff(sdiag) ; [~,L] = max(sdd) ;  %[dummy,L] = max(sdd) ; % plot(sdd) ;
  end
  
% calculate poles (eigenvalues)
  p = eig(inv(s(1:L,1:L))*u(:,1:L)'*A(:,1:K)*v(:,1:L)) ;