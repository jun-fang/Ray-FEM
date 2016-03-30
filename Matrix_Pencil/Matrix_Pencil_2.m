function [z,A] = Matrix_Pencil_2(y, tol, L)
% function [B, A, F] = mpencil(Y, X, [TOL, L])
% MPENCIL Fits a sum of exponential functions to data by a Matrix Pencil method
%
%	B = MPENCIL(Y,X) approximates Y(i) by a weighted sum of
%	exponentials sum_j=1^p A(j)*exp(B(j)*X(i)). The X(i) are
%	assumed equally spaced.
%
% The number of exponential is defined with the third optional
% parameter TOL. If TOL < 1 (and > 0) the number of terms is the same
% as the number of such singular values (of the data) that are larger
% than TOL*biggest singular value. If TOL >= 1 and an integer, the
% number of exponential terms is TOL. The default TOL value is 1e-3.
%
% The last optional parameter L is a parameter for the Matrix Pencil
% method. Values between N/3 and N/2 are preferred where N is the
% number of data (Y).
% The default L is the ceil(1/2*(ceil(N/3)+floor(N/2))).
%
%	The A(i) and the fitted values F(i) are optionally returned by
% [B,A,F] = MPENCIL(Y,X).

%	References:
% Sarkar, T.K. and Pereira, O.
% Using the Matrix Pencil Method to Estimate the Parameters of a Sum
% of Complex Exponentials.
% IEEE Antennas and Propagation Magazine, Vol. 37, No 1, Feb. 1995

%	Matti Taskinen, Rolf Nevanlinna Institute, University of Helsinki
% 06.03.1996. Latest revision 01.04.1996.


% Check the vector sizes: 
N = size(y, 1);

if size(y, 2) ~= 1
  error('Y should be column vector');
end


% Fill in the missing parameters:
if nargin < 3
  L1 = ceil(1/3 * N);
  L2 = floor(1/2 * N);
  L = ceil((L1 + L2) / 2);
  if nargin < 2
    tol = 1e-3;
  end
end

% Check the parameter values:
M_given = round(tol) == tol;
if (tol < 0) || ((~M_given) && (tol > 1))
  error('TOL should be either >= 1 and an integer or < 1 and > 0');
end

if L > N/2
  error('L shoud be < N/2');
end

if M_given && (L < tol)
  error('TOL should be <= L');
end


%% Assembling Y, nice job!
Y = zeros(N-L, L+1);
ind = 0:N-L-1;
for j = 1:L+1
  Y(:,j) = y(ind+j);
end

        
[U,S,V] = svd(Y,0); % Economy size!

if M_given
  % The number of exponential terms may be given in TOL:
  M = tol;
else
  % Otherwise figure it out from the singular values:
  D = diag(S);
  for M = 1:length(D)-1
    if (abs(D(M+1)/D(1)) <= tol)
      break;
    end
  end
end

SM = S(:,1:M);

VM = V(:,1:M);
V1 = VM(1:L,:);
V2 = VM(2:L+1,:);

Y1 = U*SM*V1';
Y2 = U*SM*V2';

A = pinv(Y1)*Y2;
z = sort(eig(A), 'descend');
% z = eig(A);
z = z(1:M);

    
