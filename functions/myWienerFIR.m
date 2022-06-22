function [w, E] = myWienerFIR(x,d,p)
%MYWIENERFIR create a non-adaptive noise-removal FIR filter
%   w = MYWIENERFIR(x,d,p) compute coefficients w of the filter based
%   on the signals x and d.
%   x and d have to be jointly wide-sense stationary (WSS).
%   x : noisy signal measured
%   d : signal that we would like to have (without noise)
%   p : filter order
%   
%   [w,E] = MYWIENERFIR(x,d,p) compute coefficients w and return the
%   estimated minimum mean-square error E.
%  
% A.Rey (c) MSE 2022 

if nargin~=3
    error('Please provide 3 arguments to the myWienerFIR function.');
end

% autocorrelation of x
rxx = xcorr(x, x, p, 'biased');
rxx = rxx(p+1:end);     % take only positive part

% cross-correlation between d and x
rdx = xcorr(d, x, p, 'biased');
% take only positive part and make it vertical
rdx = rdx(p+1:end);
rdx = rdx(:);

rdd0 = dot(d,d) / numel(d);

Rxx = toeplitz(rxx, rxx');
w = Rxx \ rdx; % inv(Rxx) * rdx;

E = rdd0 - rdx' * w(:);
end