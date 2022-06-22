function J = computeJFigure(x,a1,a2)
% function to compute the average squared error J for prediction based on the coefficients a
% x[n] = -a1*x[n-1]-a2*x[n-2]

% for simplicity, we will avoid the border of the vector x: xpred will be
% predict from k=[3 to end]
% xpred = -a1*x[2:end-1] - a2*x[1:end-2]
% e = x[3:end] - xpred

x = x(:);

% vector of x predicted
xpred = a1*x(2:end-1) + a2*x(1:end-2);
% error for each x
e = x(3:end) - xpred;

J = mean(e.^2, 1);
end