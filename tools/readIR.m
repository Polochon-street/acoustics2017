function [h,fs] = readIR(fname,fs)
%readIR   Read and resample an impulse response stored as a wave file
% 
%USAGE
%   [h,fs] = readIR(fname)
%   [h,fs] = readIR(fname,fs)
% 
%INPUT ARGUMENTS
%   fName : impulse response file name 
%      fs : requested sampling frequency in Hertz, the impulse response
%           will be resampled accordingly (default, fs = [])
% 
%OUTPUT ARGUMENTS
%    h : impulse response [nSamples x nChannels]
%   fs : sampling frequency in Hertz
% 
%   readIR(...) plots the impulse response in a new figure.
%
%   See also readAudio.

%   Developed with Matlab 9.3.0.713579 (R2017b). Please send bug reports to
%   
%   Author  :  Tobias May, © 2017
%              Technical University of Denmark
%              tobmay@elektro.dtu.dk
%
%   History :
%   v.0.1   2017/09/08
%   ***********************************************************************


%% CHECK INPUT ARGUMENTS  
% 
% 
% Check for proper input arguments
if nargin < 1 || nargin > 2
    help(mfilename);
    error('Wrong number of input arguments!')
end

% Set default values
if nargin < 2 || isempty(fs); fs = []; end


%% LOAD IMPULSE RESPONSE
% 
% 
% Check if sampling frequency was specified
if isempty(fs)
    % No sampling frequency specified, return original fs
    [h,fs] = readAudio(fname);
else
    % Resample impulse response using specified sampling frequency
    h = readAudio(fname,fs);
end


%% PLOT IMPULSE RESPONSE
%
%
% If no output is specified
if nargout == 0
    figure;
    plot(0:size(h,1)-1,h,'color',[0 0.3895 0.9712])
    xlabel('$n$','interpreter','latex')
    ylabel('$h[n]$','interpreter','latex')
    grid on;
    ylim([-1 1])
    xlim([0 size(h,1)-1])
end
