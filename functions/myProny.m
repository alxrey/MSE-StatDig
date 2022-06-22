function [ap,bq,xhat,Els] = myProny(x,p,q)
%MYPRONY computes the Prony approximation
%   [ap,bq,xhat,Els] = myProny(x,p,q) computes the Prony approximation of
%   the coefficients of an ARMA system with P poles and Q zeros, applied on
%   the signal X.
%   The system function is of the form
%          H(z) = B(z)/A(z)
% Outputs:
%   AP   : Denominator coefficients.
%   BQ   : Numerator coefficients.
%   XHAT : Approximated signal.
%   ELS  : L2 norm of the approx error x-xhat.
%
% A.Rey (c) MSE 2022

narginchk(3,3);

% force x to be in column
x = x(:);
N = length(x);

% create convolution matrix
X = convmtx(x, p+1);

% extract Xq and xq_1
Xq   = X(q+1:N+p-1, 1:p);
xq_1 = X(q+2:N+p, 1);

ap = -pinv(Xq)*xq_1;
ap = [1; ap];

X0 = X(1:q+1, 1:p+1);
bq = X0*ap;

xhat = impz(bq, ap, length(x));
Els = sqrt(sum((x-xhat).^2));

end