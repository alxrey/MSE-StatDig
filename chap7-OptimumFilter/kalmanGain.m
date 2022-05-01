function [K,xHat,P] = kalmanGain(F, H, Qw, Qv, xHatInit, PInit, z)
%MYKALMANGAIN
% This function implements the stationnary Kalman filter with
% State equation : x[n]=Fx[n-1]+w[n]
% Measurement equation : z[n]=Hx[n]+v[n]
% F is the state transition matrix
% H is the measurement matrix
% Qw is the covariance matrix of w[n]
% Qv is the covariance matrix of v[n]
% xHatInit init state vector
% PInit init covariance error
% z measurement vector z[n]

N = length(z);

%B=0;
%u=zeros(1,N);

[q, p] = size(H);
% Nobs = q
% Nstates = p

P = zeros(p, p, N);
P(:,:,1) = PInit;
K = zeros(p, q, N);
xHat = zeros(p, N);
xHat(:,1) = xHatInit;
I = eye(p);

for k=2:N
    % Prediction
    xHatk1 = F*xHat(:,k-1); % + B*u(n);
    Pk1 = F*P(:,:,k-1)*F' + Qw;

    % Update
    %K(:,:,k) = Pk1*H'*inv(H*Pk1*H' + Qv);
    % to avoid warning it's better to use the below one:
    K(:,:,k) = Pk1*H'/(H*Pk1*H' + Qv);
    P(:,:,k) = (I - K(:,:,k)*H) * Pk1;
    xHat(:,k) = xHatk1 + K(:,:,k) * (z(:,k)-H*xHatk1);
end
end