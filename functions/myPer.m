function Px = myPer(x,NFFT)
%MYPER Estimates the spectrum of a process using the periodigram.
%   PX = MYPER(X,NFFT)
%Inputs:
%   X    : input sequence
%   NFFT : number of points used for the DTFT
%Outputs:
%   PX   : estimated power spectal density using a linear scale
%
% see also myPER.
%
% A.Rey (c) MSE 2022

narginchk(1,2);

x = x(:);
N = length(x);

if nargin < 2
    FFT = fft(x);  
else
    FFT = fft(x, NFFT);
end

Px  = 1/N * abs(FFT(1:NFFT/2)).^2;

end