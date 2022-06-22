function [ap,bq,Els] = myShanks(x,p,q)
%MYSHANKS computes the Shanks approximation of a filter of ARMA system with
% P poles and Q zeros, applied on the signal HN.
%   
% Outputs:
%   AP   : Denominator coefficients.
%   BQ   : Numerator coefficients.
%   ELS  : L2 norm of the approx error x-xhat.
%
% A.Rey (c) MSE 2022

narginchk(3,3)

% force x to be in column
x = x(:);
N = length(x);

ap = myProny(x,p,q);
u = [1; zeros(N-1,1)];  % unit sample
g = filter(1, ap, u);
G = convm(g,q+1);
bq = G(1:N,:)\x;
err = x'*x-x'*G(1:N,1:q+1)*bq;
Els = norm(err);

end