function [W, E] = myRLS(x,d,nord,lambda)
% MYRLS	Adaptive filtering using the Recursive Least Squares algorithm.
%   [W,E] = myRLS(x,d,nord,lambda)
% Inputs:
%       x      : input data to the adaptive filter.
%       d      : desired output.
%       nord   : number of filter coefficients.
%       lambda : exponential forgetting factor.
% Outputs:
%       W      : output matrix of the filter coefficients.
%                n'th row contains the filter coefficients at time n.
%                m'th column contains the m'th filter coeff vs. time.
%       E      : output vector of the error sequence vs time
%
% Seel also MYLMS, MYNLMS.
%
% A.Rey (c) MSE 2022

% Inspired from:
% https://www.mathworks.com/matlabcentral/fileexchange/2183-statistical-digital-signal-processing-and-modeling

if nargin < 3
    error 'Not enough argument provided'
elseif nargin < 4
    % default value
    lambda = 1;
end

% Convolution matrix of x
X = convmtx(x, nord);

% Initialize output matrix where N is the number of samples
%N = size(X,1);
N = length(x);
W = zeros(N, nord);
E = zeros(N, 1);
% initial error
E(1) = d(1);

% Set the initial filter coefficient to 0
W(1,:) = zeros(1,nord);

% Inverse of the autocorrelation matrix
% -> a solution is to initialize the autocorrelation matrix with a simple
% constant Rx(0) = 0.001*Id
P = eye(nord)/0.001;

for k=2:N
    z = P*X(k,:)';
    g = z / (lambda + X(k,:)*z);
    alpha = d(k) - X(k,:)*W(k-1,:)';
    W(k,:) = W(k-1,:) + alpha*g';
    P = (P - g*z') / lambda;
    E(k) = d(k) - W(k-1,:)*X(k,:)';   % e[n] = d[n] - dhat[n]
end

end