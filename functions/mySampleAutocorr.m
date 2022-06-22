function [rxhat,k] = mySampleAutocorr(xn, scaleopt)
%MYSAMPLEAUTOCORR computes the autocorrelation of XN
%   [RXHAT,K] = MYSAMPLEAUTOCORR(XN), where XN is a length N vector,
%   returns the length 2*N-1 auto-correlation sequence RXHAT and a vector
%   of lag indices (K).
%
% A.Rey (c) MSE 2022

narginchk(1,2);

xn = xn(:)';
N = length(xn); % input signal length
K = 2*N - 1;    % output signal length

% preallocation
rxhat = zeros(1, K);
k = -(N-1):N-1;

xn = [xn zeros(1, K-N)]; % 0-padding

for i = 1:K
    rxhat(i) = sum(xn .* circshift(xn, k(i)));
end

% normalization if asked
if strcmp(scaleopt, 'biased')
    rxhat = rxhat ./ N;
end

end