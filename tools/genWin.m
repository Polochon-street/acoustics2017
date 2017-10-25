function w = genWin(N,type,opt)
% Create generalized cosine windows
%
%USAGE
%   w = genWin(N)
%   w = genWin(N,type,opt)
%
%INPUT ARGUMENTS
%      N : window size in samples
%   type : string specifying window type (default, type = 'hamming')
%          'rect'      = rectangular window
%          'hann'      = hann window
%          'hamming'   = hamming window 
%          'blackmann' = blackmann window
%    opt : string specifying window property (default, opt = 'periodic')
%          'periodic'  = unique maximum at N/2+1
%          'symmetric' = symmetric start and end-points
%
%OUTPUT ARGUMENTS
%   w : window function [N x 1]
% 
%   genWin(...) plots the window function in a new figure.
%
%   See also cola, iscola, stft and istft.

%   Developed with Matlab 9.2.0.538062 (R2017a). Please send bug reports to
%
%   Author  :  Tobias May, © 2017
%              Technical University of Denmark
%              tobmay@elektro.dtu.dk
%
%   History :
%   v.0.1   2017/09/07
%   ***********************************************************************


%% CHECK INPUT ARGUMENTS
%
%
% Check for proper input arguments
if nargin < 1 || nargin > 3
    help(mfilename);
    error('Wrong number of input arguments!')
end

% Set default values
if nargin < 3 || isempty(opt);  opt  = 'periodic'; end
if nargin < 2 || isempty(type); type = 'hamming';  end
    

%% CREATE GENERALIZED COSINE WINDOW
% 
% 
% Check if window should be periodic or symmetric
switch(lower(opt))
    case 'periodic'
        % Create periodic window
        N = N + 1;
        
        % Set flag
        bPeriodic = true;
    case 'symmetric'
        % Set flag
        bPeriodic = false;
    otherwise
        error('Window option "%s" is not supported.',lower(opt));
end

% Select window type
switch lower(type)
    case 'rect'
        a = 1; 
        b = 0; 
        c = 0;
    case 'hann'
        a = 0.5; 
        b = 0.5; 
        c = 0;
    case 'hamming'
        a = 0.54; 
        b = 0.46; 
        c = 0;
    case 'blackmann'
        a = 0.42; 
        b = 0.5; 
        c = 0.08;
    otherwise
        error('Window type "%s" is not supported',lower(type));
end

% Sample indices
n = (0:N-1)';

% Generalized form of cosine windows
w = a - b * cos(2*pi*n/(N-1)) + c * cos(4*pi*n/(N-1));

% For periodic windows ...
if bPeriodic
    % Delete last element (since we added N + 1 at the beginning)
    w(end) = [];
    n(end) = [];
end


%% PLOT WINDOW FUNCTION
%
%
% If no output is specified
if nargout == 0
    figure;
    hold on;
    plot(n,w,'color',[0 0.3895 0.9712])
    legend(lower(type))
    xlabel('$n$','interpreter','latex')
    ylabel('$w[n]$','interpreter','latex')
    grid on;
    xlim([n(1) n(end)] + 1)
end
