clear
close all
clc

% Octave-specific stuff
% pkg load signal

% Install subfolders
addpath irs
addpath signals
addpath tools

%% USER PARAMETERS
% 
% 
% Sampling frequency
fsHz = 16E3;

% Source signal
fileName = 'signals/speech@24kHz.wav';

% Impulse response
roomName = 'irs/church_0.wav';
% roomName = 'church_30.wav';
% roomName = 'church_270.wav';
% roomName = 'Cortex_45deg.wav';
% roomName = 'Room_A_45deg.wav';
% roomName = 'Room_B_45deg.wav';
% roomName = 'Room_C_45deg.wav';
% roomName = 'Room_D_45deg.wav';

%% CREATE BINAURAL SIGNAL
% 
% 
% Load source signal
s = readAudio(fileName,fsHz);

% Zero-padd signal
s = cat(1,s,zeros(fsHz,1));

% Load impulse response
h = readIR(roomName,fsHz);
hL = h(:, 1);
hR = h(:, 2);
% Fast convolution using the overlap-save method
xL = convolveFFT_OLS(s, hL);
xR = convolveFFT_OLS(s, hR);

tauSec = 1;
winSec = 8e-3;
fs = fsHz;

% First computations to do DDR stuff
brir = cat(2, hR, hL);
hDirect = splitBRIR(brir, fs);
hDirectL = hDirect(:, 1);
hDirectR = hDirect(:,2);
sLDry = convolveFFT_OLS(s, hDirectL);
sRDry = convolveFFT_OLS(s, hDirectR);

%% PERFORM COHERENCE-BASED DEREVERBERATION
[sL, sR, G] = dereverb(xL, xR, fs, winSec, tauSec);

%% Finding the best "tau"
DeltaDDRs = [];
tauSecs = [];
for tauSec = 10e-3:0.01:1
    [sL, sR, G] = dereverb(xL, xR, fs, winSec, tauSec);
    DDR_pre = 10*log10((sum(sLDry.*sLDry) + sum(sRDry.*sRDry)) / ...
        (sum((xL-sLDry).*(xL-sLDry) + (xR - sRDry).*(xR - sRDry))));

    DDR_post = 10*log10((sum(sLDry.*sLDry) + sum(sRDry.*sRDry)) / ...
        (sum((sL - sLDry).*(sL - sLDry) + (sR -sRDry).*(sR - sRDry))));
    DeltaDDR = DDR_post - DDR_pre;
    tauSecs = [tauSec, tauSecs];
    DeltaDDRs = [DeltaDDR, DeltaDDRs];
end

[DeltaDDR, index] = max(DeltaDDRs);
tauSec = tauSecs(index);
fprintf('Highest DeltaDDR is %f for tauSec = %f\n', DeltaDDR, tauSec);

%% PERFORMING EVALUATION

N = floor(winSec * fs); % could be floor or ceil
R = floor(N * 0.25); % could be floor or ceil
M = pow2(nextpow2(N));
w = cola(N, R, 'hamming', 'wola');

[S,t,f] = stft(s, fs, w, R, M);
[XR] = stft(xR, fs, w, R, M);
[SR] = stft(sR, fs, w, R, M);

% positive frequencies
f = f((size(f)/2)+1:end);
S = S(size(f):end,:);
XR = XR(size(f):end,:);
SR = SR(size(f):end,:);

% plots
figure;
imagesc(t,f,10*log10((abs(S).^2)/(abs(max(max(S)))^2)));
caxis([-60,0]);
colorbar;
title('Dry');
xlabel('time, s'); ylabel('frequency, Hz');
figure;
imagesc(t,f,10*log10((abs(XR).^2)/(abs(max(max(XR)))^2)));
caxis([-60,0]);
colorbar;
title('Reverb');
xlabel('time, s'); ylabel('frequency, Hz');
figure;
imagesc(t,f,10*log10((abs(SR).^2)/(abs(max(max(SR)))^2)));
caxis([-60,0]);
colorbar;
title('Dereverb');
xlabel('time, s'); ylabel('frequency, Hz');

%% SIGNAL PLOT

figure;
plot(xR,'r');
hold;
plot(sR,'g');
title('Signals');
legend('Reverb','Enhanced');
f = f((size(f)/2)+1:end);

function C = estCohere(XL, XR, alpha)
    phiLL = XL .* conj(XL);
    phiRR = XR .* conj(XR);
    phiLR = XL .* conj(XR);

    b = [1 - alpha];
    a = [1, -alpha];
    phiLL = filter(b, a, phiLL);
    phiRR = filter(b, a, phiRR);
    phiLR = filter(b, a, phiLR);

    C = phiLR ./ sqrt(phiLL .* phiRR);
end

function [sL, sR, G] = dereverb(xL, xR, fs, winSec, tauSec)
    % make winSec and tauSec optional params
    N = floor(winSec * fs); % could be floor or ceil
    R = floor(N * 0.25); % could be floor or ceil
    M = pow2(nextpow2(N));

    % Window to meet COLA constraint 3.3.1
    w = cola(N, R, 'hamming', 'wola');

    % Perform stft (3.1.1)
    [XL, tSec, fs_stft] = stft(xL, fs, w, R, M);
    XR = stft(xR, fs, w, R, M);
    
    alpha = exp(-R/(tauSec * fs));
    % Estimate short-term IC (3.4.1)
    C = estCohere(XL, XR, alpha);

    % Derive/Apply gain function
    G = abs(C) .* abs(C);

    SL = abs(XL) .* G .* exp(j * angle(XL));
    SR = abs(XR) .* G .* exp(j * angle(XR));

    % Apply ISTFT to reconstruct signal
    sL = istft(SL, w, R, 'wola', numel(xL));
    sR = istft(SR, w, R, 'wola', numel(xR));
end
