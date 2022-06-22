function [ap,bq,xhat,Els] = myPade(x,p,q)
%MYPADE computes the Pade approximation
%   [AP,BQ,ELS,XHAT] = MYPADE(X,P,Q), computes the Pade approximation of
%   the coefficients of an ARMA filter with P poles and Q zeros, applied on
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

% create convolution matrix
X = convmtx(x, p+1);

% extract Xq and xq_1
Xq   = X(q+2:q+p+1, 2:p+1);
xq_1 = X(q+2:q+p+1, 1);

if det(Xq)
    ap = Xq\-xq_1;
    ap = [1; ap];
else
   error 'Xq is a singular matrix !'
end

X0 = X(1:q+1, 1:p+1);
bq = X0*ap;
%bq = [bq; zeros(length(ap)-length(bq), 1)]; % pad bq with 0

xhat = impz(bq, ap, length(x));
%Els = sqrt(sum((x-xhat).^2));
Els = norm(x-xhat);

end