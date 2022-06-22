function [ap,bq0,ep] = covm(x,p)
%COVM Computes the all-pole model using the covariance method.
%   [AP,BQ0,ERR] = COVM(X,P) computes an all-pole model with P poles based 
%   on signal X. The model is of the form H(z) = b(0)/A(z).
% Outputs:
%     AP  : P denominator coefficients a=[1, a(1), ..., a(p-1)] 
%     BQ0 : unique coefficient b(0).
%     EP  : mean squared error
%
%
%  A.Rey (c) MSE 2022

% source:
%   Course MA-StatDig, chapter 4 (Signal approximation), page 25

narginchk(2,2);

x = x(:);
N = length(x);

if p >= N, error 'Model order to large', end

% eq. 4.130
X = convmtx(x, p+1);
Xp = X(p:N-1, 1:p);
xp = X(p+1:N, 1);

ap = [1; -pinv(Xp)*xp];
bq0 = 1;

% eq 4.126
en = conv(x, ap);
ep = sum(en(p+1:N).^2); % epsilon_p

end