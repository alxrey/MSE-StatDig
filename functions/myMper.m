function Px = myMper(x,NFFT,window)
%MYMPER	Estimates the spectrum of a process using the modified periodigram.
%	PX = MYMPER(X,NFFT,WINDOW) returns the spectrum using by default a
%   rectangular window.
%Inputs:
%   X      : input sequence
%   NFFT   : number of points used for the DTFT
%   window : indicates to use a kind of window which gives the
%            possibility to smooth the periodiagram (optional)
%Outputs
%   PX   : estimated power spectal density using a linear scale.
%
% see also myPER.
%
% A.Rey (c) MSE 2022

narginchk(1,3);

x = x(:);
N = length(x);

w  = ones(N,1);

if strcmp(window, 'hamming')
    w  = hamming(N);
elseif strcmp(window, 'bartlett')
    w = bartlett(N);
end

if nargin < 2
    FFT = fft(x.*w);
else
    FFT = fft(x.*w, NFFT);
end
Px  = 1/N * abs(FFT(1:NFFT/2)).^2;

end