function w = cola(N,R,winType,method,bPlot)
%cola   Normalized window to enforce constant overlap add criterion
%
%USAGE
%   w = cola(N)
%   w = cola(N,R,winType,method)
%
%INPUT ARGUMENTS
%         N : window size in samples
%         R : shift between adjacent windows in samples 
%             (default, R = round(numel(w)/4))
%   winType : string specifying window type (default, type = 'hamming')
%             'rect'      = rectangular window
%             'hann'      = hann window
%             'hamming'   = hamming window 
%             'blackmann' = blackmann window
%    method : string specifying reconstruction method
%             'ola'  = overlap-add method: w is applied in the analysis
%                      stage, reconstruction is based on overlapping and
%                      adding individual time frames  
%             'wola' = weighted overlap-add method: w is applied both in
%                      the analysis and synthesis stage, reconstruction is
%                      based on windowing (weighting), overlapping and
%                      adding individual time frames   
%     bPlot : plot individual and accumulated analysis/synthesis windows
%             (default, bPlot = false)
%
%OUTPUT ARGUMENTS
%   w : window function [N x 1]
% 
%   cola(...) plots the COLA criterion in a new figure.
%
%   See also iscola, genWin, stft and istft.

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
if nargin < 1 || nargin > 5
    help(mfilename);
    error('Wrong number of input arguments!')
end

% Set default values
if nargin < 2 || isempty(R);       R       = round(N/4); end
if nargin < 3 || isempty(winType); winType = 'hamming';  end
if nargin < 4 || isempty(method);  method  = 'wola';     end
if nargin < 5 || isempty(bPlot);   bPlot   = false;      end

% Check if N and R are integers
if rem(N,1) || rem(R,1)
    error('Window size "N" and step size "R" must be integers.')
end


%% CREATE WINDOW FUNCTION
% 
% 
% Create periodic window function
w = genWin(N,winType,'periodic');

% Select method
switch(lower(method))
    case 'ola'
        % Normalization constant
        K = R / sum(w);
        
    case 'wola'
        % Take square-root of window function 
        w = sqrt(w);
        
        % Normalization constant
        K = sqrt(R / sum(w.^2));
    otherwise
        error('Reconstruction method "%s" is not supported.',method);
end

% Normalized window
w = w * K;

% Check COLA criterion
iscola(w,R,method);


%% PLOT WINDOW FUNCTION
%
%
% If no output is specified
if nargout == 0 || bPlot == true
    
    % Plot COLA criterion
    iscola(w,R,method,true);
end
