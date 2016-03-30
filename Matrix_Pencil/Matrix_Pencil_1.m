function [z,A] = Matrix_Pencil_1(y, L)
%% Matrix Pencil Method
%
% Reference: Matrix Pencil Method for Estimating Parameters of bxponentially DampedKJndamped Sinusoids in Noise
% Yingo Hua and Tanpan K. Sarkar.

% Check the vector sizes: 
N = size(y, 1);

if size(y, 2) ~= 1
  error('Y should be column vector');
end


% Fill in the missing parameters:
if nargin < 2
  L1 = ceil(1/3 * N);
  L2 = floor(1/2 * N);
  L = ceil((L1 + L2) / 2);
end

if L > N/2
  error('L shoud be < N/2');
end


%% Assembling Y, nice job!
Y = zeros(N-L, L+1);
ind = N-L-1:-1:0;
for j = 1:L+1
  Y(:,j) = y(N+1 - (ind+j));
end

Y0 = Y(:,2:end);
Y1 = Y(:,1:end-1);
        
A = pinv(Y0)*Y1;
z = sort(eig(A), 'descend');

