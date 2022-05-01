function [ap,bq,Els,xhat] = pade(x,p,q)
%PADE compute the pade approximation of the coeff of ARMA system with p 
%poles and q zeros usage: [ap,bq,Els,xhat] = PADE(x,p,q)
%Els is optional and is the L2 norm of the approx error x-xhat
%xhat is optional and it is the approximated signal.
%    (c) A.Rey MSE 2022 r1.0

x = x(:); % force x to be in column

X = convmtx(x, p+1);

Xq = X(q+2:q+p+1, 2:p+1);
xq_1 = X(q+2:q+p+1, 1);

if det(Xq)
    ap = inv(Xq)*(-xq_1); % pinv is used in case
    ap = [1; ap];
else
    disp('Xq is a singular matrix !')
    return
end

X0 = X(1:q+1, 1:p+1);
bq = X0*ap;

xhat = impz(bq, ap, length(x));
Els = norm(x-xhat);

end