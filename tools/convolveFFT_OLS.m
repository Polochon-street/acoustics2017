function y = convolveFFT_OLS(x,h,N)
%convolveFFT_OLS   Overlap-save convolution in the frequency domain.
% 
%USAGE
%   y = convolveFFT_OLS(x,h)
%   y = convolveFFT_OLS(x,h,N)
% 
%INPUT ARGUMENTS
%   x : input sequence [Nx x 1]
%   h : impulse response [M x 1]
%   N : FFT size (default, N = 2 * pow2(nextpow2(M)))
% 
%OUTPUT ARGUMENTS
%   y : output sequence [Nx x 1]
% 
%   See also convolve, convolveFFT and convolveFFT_OLA.

%   Developed with Matlab 9.2.0.538062 (R2017a). Please send bug reports to
%   
%   Author  :  Tobias May, © 2017
%              Technical University of Denmark (DTU)
%              tobmay@elektro.dtu.dk
%
%   History :
%   v.0.1   2017/08/12
%   ***********************************************************************


%% CHECK INPUT ARGUMENTS
%
%
% Check for proper input arguments
if nargin < 2 || nargin > 3
    help(mfilename);
    error('Wrong number of input arguments!')
end

% Get dimensions
xDim = size(x);
hDim = size(h);

% Check if x and h are vectors
if min(xDim) > 1 || min(hDim) > 1
    error('x and h must be vectors.')
else
    % Ensure column vectors
    x = x(:);
    h = h(:);
    
    % Dimensionality of x and h
    Nx = max(xDim);
    M = max(hDim);
end

% Set default values
if nargin < 3 || isempty(N); N = 2^(nextpow2(M) + 1); end


%% OBTAIN SPECTRUM OF FIR FILTER
% 
% 
% Block size
L = N - M + 1;

% Zero-pad impulse response
h = cat(1,h,zeros(N - M,1));

% Compute H[k]
H = fft(h);


%% BLOCK-BASED CONVOLUTION IN THE FREQUENCY DOMAIN
% 
% 
% Number of blocks
nBlocks = ceil(Nx / L);

% Zero-pad input such that an integer number of blocks can be processed
x = cat(1,zeros(M-1,1),x,zeros(nBlocks * L - Nx,1));

% Allocate memory
y = zeros(nBlocks * L,1);

% Loop over the number of blocks
for ii = 1 : nBlocks

    % Extract ii-th block (L input samples plus M - 1 preceeding samples)
    xm = x((1:L + M - 1) + (ii-1) * L);
    
    % DFT
    Xm = fft(xm);
    
    % Multiplication with impulse response spectrum
    Ym = Xm .* H;
    
    % IDFT
    ym = ifft(Ym);
        
    % Overlap-save L output samples (discard first M - 1 points)
    y((1:L) + (ii-1) * L) = ym(M:L + M - 1);
end

% Trim output signal
y = y(1:Nx);

