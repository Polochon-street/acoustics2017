function [bCola,offset] = iscola(w,R,method,bPlot)
%iscola   Check if window satisfies constant overlap add criterion.
%
%USAGE
%   cola(w,R,method)
%   [bCola,offset] = cola(w,R,method)
%   [bCola,offset] = cola(w,R,method,bPlot)
%
%INPUT ARGUMENTS
%         w : analysis/synthesis window function [N x 1]
%         R : shift between adjacent windows in samples 
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
%    bCola : binary flag indicating if COLA criterion is satisfied
%   offset : sample offset after which perfect reconstruction is achieved
% 
%   iscola(...) will display an error message if the COLA criterion is not
%   satisfied.
%
%   See also cola and genWin.

%   Developed with Matlab 9.3.0.713579 (R2017b). Please send bug reports to
%
%   Author  :  Tobias May, © 2017
%              Technical University of Denmark
%              tobmay@elektro.dtu.dk
%
%   History :
%   v.0.1   2017/10/10
%   ***********************************************************************


%% CHECK INPUT ARGUMENTS
%
%
% Check for proper input arguments
if nargin < 3 || nargin > 4
    help(mfilename);
    error('Wrong number of input arguments!')
end

% Set default values
if nargin < 4 || isempty(bPlot); bPlot = false; end

% Check if R is an integer
if rem(R,1)
    error('Window step size "R" must be an integer.')
end
  

%% PERFORM OLA OR WOLA RECONSTRUCTION
% 
% 
% Determine window length
N = numel(w);

% offset => perfect reconstruction
offset = (N/R-1) * R + 1;

% Transition, 2 windows, Transition
T = round(((N/R-1)*2+2)*R);

% Number of overlapping windows
L = floor((T - (N - R))/R);

% Allocate memory
z = zeros(T,L); 

% Loop over the number of frames
for ii = 1 : L
    
    % Window indices
    wIdx = (1:N) + (ii-1) * R;
    
    % Select reconstruction method
    switch(lower(method))
        case 'ola'
            % OLA: Analysis window only
            z(wIdx,ii) = w;
        case 'wola'
            % WOLA: Analysis and synthesis window
            z(wIdx,ii) = w.^2;
        otherwise
            error('Reconstruction method "%s" is not supported.',...
                lower(method))
    end
end


%% CHECK COLA CRITERION
% 
% 
% Total sum of all windows
y = sum(z,2);

% RMS erros
rmsError = rms(y(offset:end-offset+1)-1);

% Check if COLA constraint is satisfied
bCola = rmsError < 1E-10;


%% DISPLAY ERROR
% 
% 
% Display error message
if nargout == 0 && ~bCola
    error(['COLA criterion is not satisfied. Change the window ',...
        'length "N", the window function "w" or the step size "R".'])
end


%% PLOT WINDOWS AND COLA CRITERION
%
%
% Plot if specified
if bPlot
    figure;
    hold on;
    h1 = plot(z,'color',[0 0.3895 0.9712]);
    h2 = plot(y,'k-','linewidth',1);
    if bCola
        h3 = plot([offset offset; offset T-offset+1; ...
            offset T-offset+1; T-offset+1 T-offset+1]',...
            [0 1;0 0;1 1; 0 1]','--','color',[0 0.5 0],'linewidth',2);
        legend([h1(1) h2 h3(1)],{'$w[n]$' '$\sum\limits w[n]$' ...
            'COLA = 1'},'interpreter','latex')
        ylim([0 1])
    else
        legend([h1(1) h2],{'$w[n]$' '$\sum\limits w[n]$'},...
            'interpreter','latex')
    end
    xlim([1 T])
    xlabel('$n$','interpreter','latex')
    ylabel('Amplitude','interpreter','latex')
    grid on;
end