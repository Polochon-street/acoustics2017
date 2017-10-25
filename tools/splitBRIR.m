function [d,r] = splitBRIR(brir,fsHz,timeDirect,maxITD,thresdB,bLink)
%splitBRIR   Split a BRIR into a direct and a reverberant part.
%
%USAGE
%   [d,r] = splitBRIR(brir,fsHz)
%   [d,r] = splitBRIR(brir,fsHz,timeDirect,maxITD,thresdB)
%
%INPUT PARAMETERS
%         brir : binaural room impulse response [nPoints x 2] 
%         fsHz : sampling frequency in Hertz
%   timeDirect : duration of the direct part starting from the maximum peak
%                position in each channel (default, timeDirect = 2.5E-3)
%                According to [1], 2.5E-3 reflects the average duration of
%                anechoically measured head-related impulse responses.  
%       maxITD : maximum time delay in seconds between the main peaks 
%                in the left and right BRIR (default, maxITD = 1E-3)
%      thresdB : remove samples preceeding the maximum peak from the direct
%                part of the BRIR if their energy is below "thresdB" with
%                respect to the maximum peak value (default, thresdB = inf)
% 
%OUTPUT PARAMETERS
%   d : direct part of the BRIR [nPoints x 2] 
%   r : reverberant part of the BRIR [nPoints x 2] 
% 
%   splitBRIR(...) plots the BRIRs in a new figure. 
% 
%   See also filterBRIR. 
% 
%REFERENCES
%   [1] P. Zahorik, "Direct-to-reverberant energy ratio sensitivity", The
%       Journal of the Acoustical Society of America, 112(5), pp.2110-2117,
%       2002.   

%   Developed with Matlab 9.0.0.341360 (R2016a). Please send bug reports to
%   
%   Author  :  Tobias May, ? 2016
%              Technical University of Denmark
%              tobmay@elektro.dtu.dk
% 
%   History :  
%   v.0.1   2016/09/15
%   ***********************************************************************


%% CHECK INPUT ARGUMENTS  
% 
% 
% Check for proper input arguments
if nargin < 2 || nargin > 6
    help(mfilename);
    error('Wrong number of input arguments!')
end

% Set default values
if nargin < 3 || isempty(timeDirect);  timeDirect = 2.5E-3; end
if nargin < 4 || isempty(maxITD);      maxITD     = 1E-3;   end
if nargin < 5 || isempty(thresdB);     thresdB    = inf;    end
if nargin < 6 || isempty(bLink);       bLink      = false;  end


%% CREATE WINDOWS TO EXTRACT THE DIRECT AND THE REVERBERANT PART
% 
% 
% Find peak for each channel
[mVal,mIdx] = max(abs(brir),[],1);

% Dimensions
[nSamples,nChannels] = size(brir);

% Maximum ITD range in samples
maxDelayRange = 1:nSamples - mIdx(argmax(mVal));

% Apply ITD constraint
if isfinite(maxITD)
    maxDelayRange(round(fsHz * maxITD)+1:end) = [];
end

if bLink
    mIdx = [min(mIdx),min(mIdx)];
else
    % Ensure that delay betwen both peaks is within predefined ITD limit
    if mVal(1) > mVal(2)
        mIdx(2) = mIdx(1) + argmax(abs(brir(mIdx(1) + maxDelayRange,2)));
    else
        mIdx(1) = mIdx(2) + argmax(abs(brir(mIdx(2) + maxDelayRange,1)));
    end
end

% Allocate memory
winD = zeros(nSamples,nChannels);

% Create window for the direct part
winD(1:min(nSamples,mIdx(1) + round(fsHz * timeDirect)),1) = 1;
winD(1:min(nSamples,mIdx(2) + round(fsHz * timeDirect)),2) = 1;

% The reverberant window is the inverse of the direct window
winR = 1 - winD;

% Remove samples preceeding the maximum peak from the direct part of the
% BRIR if their energy is below "thresdB" with respect to the maximum peak
% value 
if isfinite(thresdB)
    thres      = db2pow(-abs(thresdB));
    bSet2ZeroL = brir(1:mIdx(1),1).^2 < brir(mIdx(1),1).^2 * thres;
    bSet2ZeroR = brir(1:mIdx(2),2).^2 < brir(mIdx(2),2).^2 * thres;
    
    winD(bSet2ZeroL,1) = 0;
    winD(bSet2ZeroR,2) = 0;
end


%% APPLY WINDOWS
% 
% 
% Filter direct and reverberant part
d = brir .* winD;
r = brir .* winR;


%% PLOT BRIRs
% 
% 
% Show direct and reverberant part
if nargout == 0
    
    figure;
    subplot(211);
    hold on;
    plot(d(:,1),'k');
    plot(r(:,1),'b');
    legend({'direct' 'reverberant'})
    title('left ear')
    
    xlabel('Number of samples')
    ylabel('Amplitude')
    grid on;
    xlim([1 nSamples])
    ylim([-1 1])
    
    subplot(212);
    hold on;
    plot(d(:,2),'k');
    plot(r(:,2),'b');
    legend({'direct' 'reverberant'})
    title('right ear')
    
    xlabel('Number of samples')
    ylabel('Amplitude')
    grid on;
    xlim([1 nSamples])
    ylim([-1 1])
end


%   ***********************************************************************
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
% 
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
% 
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%   ***********************************************************************