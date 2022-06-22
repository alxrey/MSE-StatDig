function [K, xHat, P] = kalmanGain(F, H, Qw, Qv, xHatInit, PInit, z, B, u)
%KALMANGAIN Computes gain of a stationnary Kalman filter
% State equation : x[n]=Fx[n-1]+w[n]
% Measurement equation : z[n]=Hx[n]+v[n]
% F         transition matrix                                         [pxp]
% H         measurement matrix
% Qw        covariance matrix of w[n]
% Qv        covariance matrix of v[n]
% xHatInit  initial condition of the state vector
% PInit     initial condition of the covariance error matrix
% z         measurements vector
% B         control input model matrix                                [pxI]
% u         known input                                               [Ix1]
%
% A.Rey (c) MSE 2022

narginchk(7,9);

N = length(z);

if nargin==7
    B=0;
    u=zeros(1,N);
end

[q, p] = size(H);

% prealloc matrices
P = zeros(p, p, N);
K = zeros(p, q, N);
xHat = zeros(p, N);
I = eye(p);

% Set initial condition
P(:,:,1) = PInit;
xHat(:,1) = xHatInit;

for k=2:N
    % Prediction
    xHatk1 = F*xHat(:,k-1) + B*u(k);    % xHatk1 = xHat[n|n-1]
    Pk1 = F*P(:,:,k-1)*F' + Qw;         % Pk1 = P[n|n-1]

    % Update
    K(:,:,k) = Pk1*H'*inv(H*Pk1*H' + Qv);
    % to avoid warning it's better to use the below one:
    %K(:,:,k) = Pk1*H'/(H*Pk1*H' + Qv);
    P(:,:,k) = (I - K(:,:,k)*H) * Pk1;
    xHat(:,k) = xHatk1 + (K(:,:,k) * (z(:,k)-H*xHatk1));
end
end