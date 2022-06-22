function [A,E] = myLMS(x,d,mu,nord)
% MYLMS	Adaptive filtering using the Widrow-Hoff LMS algorithm.
%   [A,E] = MYLMS(x,d,mu,nord)
% Inputs:
%       x    : input data to the adaptive filter.
%       d    : desired output.
%       mu   : adaptive filtering update (step-size) parameter.
%       nord : number of filter coefficients.
% Outputs:
%       A    : output matrix of the filter coefficients.
%              n'th row contains the filter coefficients at time n.
%              m'th column contains the m'th filter coeff vs. time.
%       E    : output vector of the error sequence vs time
%
% Seel also MYNLMS, MYRLS.
%
% A.Rey (c) MSE 2022

% Inspired from:
% https://www.mathworks.com/matlabcentral/fileexchange/2183-statistical-digital-signal-processing-and-modeling

% ensure to get signals vertical
x = x(:);
d = d(:);

% Create the convolution matrix for x
% all the row will be used except the nord last ones
X = convmtx(x, nord);

M = length(x);

% Preallocate output matrix
A = zeros(M, nord);
E = zeros(M, 1);

% Set the initial filter coefficients to 0
%a0 = zeros(1, nord);

% compute first error
E(1) = d(1);
%E(1) = d(1) - a0*X(1,:)'; % e[1] = d[1] - dhat[1]

% compute first coefficient
A(1,:) = mu*E(1)*conj(X(1,:));
%A(1,:) = a0 + mu*E(1)*conj(X(1,:));

% Iter to compute the following estimated filter coefficients based on the
% LMS algorithm
if M > 1
    for k=2:M
        E(k)   = d(k) - A(k-1,:)*X(k,:)';         % e[n] = d[n] - dhat[n]
        A(k,:) = A(k-1,:) + mu*E(k)*conj(X(k,:)); % w[n+1] = w[n] + mu*e[n]*conj(x[n])
    end
end
end