function [ y ] = convolveFFT_OLS( x, h)
%CONVOLVEFFT_OLS performs the convolution of signal x
%   with impulse response h, using the overlap-save method
%   N = FFT size

M = length(h);
Nx = length(x);
N = 2*pow2(nextpow2(M));
L = N - M + 1;


H = fft(h,N);
A = ceil(Nx/L);
x = [x; zeros(A*L - Nx,1)];
x = [zeros(M-1,1);x];
y = [];

for n = 0:A-1
    
    xm = x(1+n*L:L+n*L+M-1);
    Xm = fft(xm);
    Ym = Xm.*H;
    ym = ifft(Ym);
    y = [y; ym(M-1+1:end)];
    
end
        
y = y(1:Nx);

end

