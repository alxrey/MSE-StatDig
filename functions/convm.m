function X = convm(x,p)
% CONVM  Generates a convolution matrix
%   X = CONVM(X,P) returns the convolution matrix, X, such that
%   convolution of x, of length N, and another vector, y, with the same
%   vector orientation as x may be expressed by matrix multiplication:
%   When x and y are column vectors, X*y is the same as CONV(x,y);
%   when x and y are row vectors, y*X is the same as CONV(x,y).
%
%   X is of size N+p-1 by p and has the following form
%
%              |  x(0)  0      0     ...      0    |
%              |  x(1) x(0)    0     ...      0    |
%              |  x(2) x(1)   x(0)   ...      0    |
%         X =  |   .    .      .              .    |
%              |   .    .      .              .    |
%              |   .    .      .              .    |
%              |  x(N) x(N-1) x(N-2) ...  x(N-p+1) |
%              |   0   x(N)   x(N-1) ...  x(N-p+2) |
%              |   .    .      .              .    |
%              |   .    .      .              .    |
%              |   0    0      0     ...    x(N)   |
%
%   CONVM and convmtx are the same but written myself as exercise
%   See also CONV, CONVMTX.
%
% A.Rey (c) MSE 2022

X = zeros(length(x)+p-1, p);

xpad = [x(:) ; zeros(p-1, 1)];

for i=1:p
    X(:,i) = circshift(xpad, i-1);
end

end