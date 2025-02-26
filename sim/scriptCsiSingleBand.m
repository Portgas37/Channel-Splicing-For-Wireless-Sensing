%scriptCsiSingleBand.m
%
% This script demonstrates how to estimate the channel frequency and impulse response from a single Wi-Fi frame.
% The Wi-Fi frame is passed through a multipath channel for which we have specified a number
% of multi-path components (MPCs) allowing us to compare the estimated CIR with the true one.
%
% In this script, we also consider different channel models and estimating the CSI from both the
% HT-LTF fields and the L-LTF fields.
%
% See:
%   - <https://ch.mathworks.com/help/comm/ug/fading-channels.html>
%   - <https://ch.mathworks.com/help/comm/ug/multipath-fading-channel.html>
%
% Note that we can work with different channel models
% comm.RayleighChannel | Rayleigh fading channel | One or more major reflected paths
% comm.RicianChannel   | Rician fading channel   | One direct line-of-sight path, possibly combined with one or more major reflected paths
%
% As such, a Rician channel is a good model for a line-of-sight (LOS) channel with some multipath like a monostatic radar.
% and a Rayleigh channel is a good model for a non-line-of-sight (NLOS) channel with no LOS component such as bistatic radar.

%% Settings
%------------------------------------------------------------------------------

clear; clc;
close all;
rng("default");

run('spliceSettingsCommon.m');

% Parameters
scenarioType  = 'test';     % One of ('test' | 'real')
scenarioIndex = 6;          % Multipath scenario
noiseEnable   = true;      % Add AWGN
channelType   = 'rayleigh'; % One of ('rayleigh' | 'rician' | 'none')

snrdB = 10;  % Signal-to-noise ratio in dB

symOffset = 1.0; % Symbol sampling offset as fraction of cyclic prefix length (0 to 1, inclusive).
                 % Default in MATLAB is 0.75. If you pick symOffset < 1.0 it will shift the estimated
                 % channel impulse response and potentially spread it out depending on whether
                 % symOffset is an integer multiple of the number of cyclic prefix samples or not (there are 16)

smoothingSpan = 1; % # sub-carriers to average together when estimating the CSI (assumes CFR is similar for neighouring sub-carriers.)

%% Channel
%------------------------------------------------------------------------------

[scenario, ~] = getMultipathScenario(scenarioType, scenarioIndex);
if strcmpi(channelType, 'rayleigh') || strcmpi(channelType, 'rician')
  channelObj = generateChannel(channelType, scenario, fs);
end

%% Transmit data
%------------------------------------------------------------------------------

txData = wlanWaveformGenerator([], cfgHT);

if strcmpi(channelType, 'rayleigh') || strcmpi(channelType, 'rician')
  txCh = channelObj(txData);
else
  txCh = txData;
end

if noiseEnable
  txCh = awgn(txCh, snrdB, 'measured');
end

%% Extract CSI
%------------------------------------------------------------------------------

% Extract L-LTF and HTLTF parts of the frame (training symbol) and perform channel estimation
txSLTF  = txCh(ind.LSTF(1) : ind.LSTF(2), :);
txLLTF  = txCh(ind.LLTF(1) : ind.LLTF(2), :);
txHTLTF = txCh(ind.HTLTF(1) : ind.HTLTF(2), :);

% Demodulate (this takes the OFDM symbol in time-domain and converts it to frequency domain)
txLLTFSymbol  = wlanLLTFDemodulate(txLLTF, cfgHT, symOffset);
txHTLTFSymbol = wlanHTLTFDemodulate(txHTLTF, cfgHT, symOffset);

% Estimate CSI
csiLLTF  = wlanLLTFChannelEstimate(txLLTFSymbol, cfgHT, smoothingSpan);
csiHTLTF = wlanHTLTFChannelEstimate(txHTLTFSymbol, cfgHT, smoothingSpan);

%% Estimate CIR
%------------------------------------------------------------------------------

% DTF Matrix for CIR estimation
prunedDFTMatLLTF  = getPrunedDFTMatrix(nfft, infoLLTF.ActiveFFTIndices);
prunedDFTMatHTLTF = getPrunedDFTMatrix(nfft, infoHTLTF.ActiveFFTIndices);

% CIR calculation
% Note that as we do not have the full number of rows we can't just use ifft/fft but we have to use
% a pruned version of the DFT matrix where we only keep the sub-carriers that exist for the L-LTF or HT-LTF fields
cirLLTF  = pinv(prunedDFTMatLLTF) * csiLLTF;
cirHTLTF = pinv(prunedDFTMatHTLTF) * csiHTLTF;

%% Plotting
%------------------------------------------------------------------------------

% Note the difference between L-LTF and HT-LTF. L-LTF uses 52 sub-carriers and HT-LTF uses 56  sub-carriers.
freqMHzLLTF  = infoLLTF.ActiveFrequencyIndices  * scSpacing / 1e6;
freqMHzHTLTF = infoHTLTF.ActiveFrequencyIndices * scSpacing / 1e6;

lenLLTF  = length(cirLLTF);
lenHTLTF = length(cirHTLTF);

timeSecLLTF  = (0:1:(lenLLTF-1)) / fs;
timeSecHTLTF = (0:1:(lenHTLTF-1)) / fs;

% Using the function for Frequency Responses
plotComplexResponse(freqMHzLLTF, csiLLTF, 'title', 'Est. CFR L-LTF', 'xLabel', 'Freq. (MHz)', ...
  'useDbForMag', false, 'unwrapPhase', true, 'markerStyle', '-o', 'useStem', false);
plotComplexResponse(timeSecLLTF, cirLLTF, 'title', 'Est.s CIR L-LTF', 'xLabel', 'Time (s)', ...
  'useDbForMag', false, 'unwrapPhase', false, 'markerStyle', '-o', 'useStem', true);

plotComplexResponse(freqMHzHTLTF, csiHTLTF, 'title', 'Est. CFR HT-LTF', 'xLabel', 'Freq. (MHz)', ...
  'useDbForMag', false, 'unwrapPhase', true, 'markerStyle', '-o', 'useStem', false);
plotComplexResponse(timeSecHTLTF, cirHTLTF, 'title', 'Est.s CIR HT-LTF', 'xLabel', 'Time (s)', ...
  'useDbForMag', false, 'unwrapPhase', false, 'markerStyle', '-o', 'useStem', true);

% True channel
visualizeChannel(channelObj, fs, true, true);
