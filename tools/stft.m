function [X,t,f] = stft(x,fs,w,R,M)
% Compute the short-time Fourier transform
%
%USAGE
%   [X,t,f] = stft(x,fs,w,R,M)
%
%INPUT ARGUMENTS
%    x : input signal [Nx x 1 | 1 x Nx] (x must be a vector)
%   fs : sampling frequency in Hertz
%    w : length-N vector containing the window function [N x 1 | 1 x N]
%    R : shift between adjacent windows in samples
%    M : DFT size, windows are zero-padded if M > N
%
%OUTPUT ARGUMENTS
%   X : complex STFT matrix [M x L]
%   t : time vector in seconds [1 x L]
%   f : frequency vector in Hertz [M x 1]

%Derive N and O
N = length(w);
O = N - R;

%Compute number of frames L
Nx = length(x);
L = ceil((Nx - O)/R);

%Zero-pad input signal x with O+LR zeros
if M > N
    w = [w; zeros(M-N, 1)];
end

x = [x; zeros(O+L*R-Nx, 1)];

%Allocate memory
xn = zeros(M,1);
X = zeros(M,L);

%Segmentation and windowing, DFT
for n = 1:L
    Rs = (n-1)*R;                           %Shift between adjacent windows
    xn = [x(1+Rs:N+Rs); zeros(M-N, 1)].*w;   %Segmentation and windowing
    X(:,n) = fft(xn);                        %DFT
end

t = (N/2+(0:L-1)*R)/fs;
f = (-fs/2:fs/M:fs/2-fs/M)';


end

