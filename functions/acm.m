function [ap,bq0,ep] = acm(x,p)
%ACM Computes the all-pole model using the autocorrelation method.
%   [AP,BQ0,ERR] = ACM(X,P) computes an all-pole model with P poles based 
%   on signal X. The model is of the form H(z) = b(0)/A(z).
% Outputs:
%	  AP  : P denominator coefficients a=[1, a(1), ..., a(p-1)]
%     BQ0 : unique coefficient b(0)
%     ERR : norm of the error
%
%   A.Rey (c) MSE 2022

narginchk(2,2);

x = x(:);
N = length(x);

X = convmtx(x, p+1);
Xp = X(1:N+p-1, 1:p);
Rx = Xp' * Xp / N;
rx = Xp' * X(2:N+p, 1) / N;
ap = [1; -Rx\rx];

ep = (Rx(1,1) + rx'*ap(2:3)) * N; % epsilon_p

bq0 = sqrt(err);

end