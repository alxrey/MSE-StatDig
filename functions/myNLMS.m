function [A,E] = myNLMS(x,d,beta,nord,epsilon)
% MYNLMS Normalized LMS adaptive filtering algorithm.
%   [A,E] = myNLMS(x,d,beta,nord,Epsilon)
% Inputs:
%       x       : input data to the adaptive filter.
%       d       : desired output
%       beta    : adaptive filtering update (step-size) parameter
%       nord    : number of filter coefficients
%       Epsilon : small value to avoid numerical issue during scaling
%                 (~1e-4)
% Outputs:
%       A    : filter coefficients output matrix
%              n'th row contains the filter coefficients at time n
%              m'th column contains the m'th filter coeff vs. time
%       E    : output vector containing error sequence vs time
%
% See also MYLMS, MYRLS.
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

E(1) = d(1);

% normalization factor (denominator) 1/(||x||+epsilon)
x_norm = X(1,:)*X(1,:)' + epsilon;

A(1,:) = beta*E(1)/x_norm*conj(X(1,:));

if M > 1
    for k=2:M
        E(k)   = d(k) - A(k-1,:)*X(k,:)';         % e[n] = d[n] - dhat[n]
        x_norm = X(k,:)*X(k,:)' + epsilon;
        A(k,:) = A(k-1,:) + beta*E(k)/x_norm*conj(X(k,:)); % w[n+1] = w[n] + mu*e[n]*conj(x[n])
    end
end
end
