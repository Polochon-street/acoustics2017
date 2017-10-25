function x = istft(X,w,R,method,Nx)
% Compute the inverse short-time Fourier transform
%
%USAGE
%   x = istft(X,w,R,method)
%   x = istft(X,w,R,method,Nx)
%
%INPUT ARGUMENTS
%        X : complex STFT matrix [M x L]
%        w : length-N vector containing the window function [N x 1 | 1 x N]
%        R : shift between adjacent windows in samples 
%   method : reconstruction method
%            'ola'  = overlap-add method: w is applied in the analysis
%                     stage, reconstruction is based on overlapping and
%                     adding individual time frames  
%            'wola' = weighted overlap-add method: w is applied both in the
%                     analysis and synthesis stage, reconstruction is based
%                     on windowing (weighting), overlapping and adding
%                     individual time frames    
%       Nx : original signal length of STFT representation X, if specified
%            the output signal x will be trimmed to [Nx x 1]
%
%OUTPUT ARGUMENTS
%   x : ouput signal [N + L * R x 1] | [Nx x 1]
%
%   istft(...) plots the reconstructed output signal in a new figure.
% 
%   See also stft, cola and genWin.

%   Developed with Matlab 9.2.0.538062 (R2017a). Please send bug reports to
%
%   Author  :  Tobias May, © 2017
%              Technical University of Denmark
%              tobmay@elektro.dtu.dk
%
%   History :
%   v.0.1   2017/09/18
%   ***********************************************************************


%% CHECK INPUT ARGUMENTS
%
%
% Check for proper input arguments
if nargin < 4 || nargin > 5
    help(mfilename);
    error('Wrong number of input arguments!')
end


%% SELECT RECONSTRUCTION METHOD
% 
% 
% Select method
switch(lower(method))
    case 'ola'
        % Overlap-add method, no synthesis window
        bWin = false;
    case 'wola'
        % Weighted overlap-add method, analsis and synthesis window
        bWin = true;
    otherwise
        error('Reconstruction method "%s" is not supported.',method);
end

% Check COLA criterion
iscola(w,R,method);


%% INVERSE SHORT-TIME FOURIER TRANSFORM
% 
% 
% Get dimensionality
[M,L] = size(X);

% Window size
N = numel(w);

% Total number of samples
T = N + L * R; 

% Allocate memory
x = zeros(T,1);

% Loop over the number of frames
for ii = 1:L
    % IDFT
    xm = real(ifft(X(:,ii),M));
    
    % Apply synthesis window
    if bWin
        xw = xm(1:N) .* w;
    else
        xw = xm(1:N);
    end
    
    % Overlap and add individual time frames
    x((1:N)+(ii-1)*R) = x((1:N)+(ii-1)*R) + xw;
end

% If specified, trim to specified length Nx
if exist('Nx','var')
    if T > Nx
        x(Nx+1:end) = [];
    else
        warning(['No trimming is performed because the reconstructed '...
            'signal "x" is shorter than "Nx".'])
    end
end


%% PLOT OUTPUT SIGNAL
%
%
% If no output is specified
if nargout == 0
    n = (0:numel(x)-1);
    
    figure;
    plot(n,x,'color',[0 0.3895 0.9712])
    grid on;
    xlim([n(1) n(end)])
    axis xy
    xlabel('n','interpreter','latex')
    ylabel('Amplitude','interpreter','latex')
end
