% scriptCsiMultiBand.m
%
% This script aims to demonstrate the process of combining Channel State Information (CSI)
% from multiple carriers to perform channel splicing. It builds on the simpler 'scriptCsiSingle.m'
%
% PERFORMANCE
%   Note that if you use rf modelType, the performance will be worse than with bb. This is because
%   we have to simulate at GHz frequencies.
%
% REFERENCES
%   - Interpolation/Decimation: Figure 11.12 in Software Defined Radio using MATLAB & Simulink and the RTL-SDR

%% Settings
%------------------------------------------------------------------------------

clear; clc;
close all;
rng("default");

run('spliceSettingsCommon.m');

% Parameters
scenarioType      = 'real'; % One of ('test' | 'real')
scenarioIndex     = 2;
noiseEnable       = true;
channelType       = 'rayleigh'; % One of ('rayleigh' | 'rician' | 'none')
randomPhaseOffset = true; % Add a random phase offset to the signal. This represents the random phase offset
                          % That occurs when we change the carrier frequency.

snrdB               = 20;
symOffset           = 1.0;
smoothingSpan       = 3; % MATLAB smoothing function
smoothingWindowSize = 3; % Additional smoothing

modelType = 'bb'; % One of ('bb' | 'rf'). If 'bb', the signal stays in baseband, if rf, the signal is mixed up to a carrier frequency and converted to a passband signal

plotAllCFR = false; % ( false | true ). If true, plot all intermediate CFRs used in the splicing process
plotMatlab = false; % ( false | true ). If true, plot the channel as MATLAB visualizes it. This is useful to compare the results.
                    % However, please read the documentation of visualizeChannel() to understand the limitations of this function.

% Channel hopping parameters
fc0    = 2.4e9 + 20e6; % First carrier frequency
fcStep = 16;          % Jump X carrier spacings. If fcStep=16 jump 5 MHz which is Wi-Fi channel spacing
%nFc    = 200; % Number of channels to hop
nFc =200;

fcList = fc0 + (0:(nFc-1)) * fcStep * scSpacing;


%% Splicing setup
%------------------------------------------------------------------------------

infoSplice = spliceInitInfo(infoHTLTF, scSpacing, fcList);
infoHTLTF.fftDcIdx = getFftDcIdx(infoHTLTF.FFTLength);

% A very basic approach to calculating the channel impulse response (CIR) is to use the inverse DFT
% where the DFT matrix is pruned to only include the active subcarriers. This is necessary because
% we do not have frequency values at the inactive subcarriers.
prunedDFTMatHTLTF    = getPrunedDFTMatrix(nfft, infoHTLTF.ActiveFFTIndices, infoHTLTF.fftDcIdx);
prunedIDFTMatHTLTF   = pinv(prunedDFTMatHTLTF);
prunedDFTMatSpliced  = getPrunedDFTMatrix(infoSplice.FFTLength, infoSplice.ActiveFFTIndices, infoSplice.fftDcIdx);
prunedIDFTMatSpliced = pinv(prunedDFTMatSpliced);

fsSplicedChannel = infoSplice.SampleRate;

%% Channel
%------------------------------------------------------------------------------

[scenario, requiredBandwidth] = getMultipathScenario(scenarioType, scenarioIndex);

if strcmpi(modelType, 'rf')
  interpFactor = ceil(4 * infoSplice.CarrierFrequenciesPassband(end) / fs); % Ensure higher than Nyquist
else
  interpFactor = ceil(fsSplicedChannel / fs);
end

fsCh = interpFactor * fs;

if strcmpi(channelType, 'rayleigh') || strcmpi(channelType, 'rician')
  channelObj = generateChannel(channelType, scenario, fsCh);
end

%% Simulate
%------------------------------------------------------------------------------


% We measure the wider channel by upconverting the narrow BW signal
txData      = wlanWaveformGenerator([], cfgHT);
txDataArray = zeros(nFc, length(txData), 'like', 1j);
rxDataArray = zeros(nFc, length(txData), 'like', 1j);

csiSegmentsArray = zeros(nFc, infoSplice.infoWLANField.NumTones, 'like', 1j); % Data to splice bundled in bursts

for iFc = 1:nFc
  txDataArray(iFc, :) = txData; % We just repeat the same signal for all channels

  if strcmpi(modelType, 'rf')
    fc = infoSplice.CarrierFrequenciesPassband(iFc);
  else
    fc = infoSplice.CarrierFrequenciesBaseband(iFc);
  end

  % Resample data for upconversion to the channel frequency
  txDataUpsampled = resample(txData, fsCh, fs);

  if strcmpi(modelType, 'rf')
    % Lowpass filter: Mostly needed when upsample/downsample is used, resample() includes a lowpass filter
    txDataUpsampledReal = lowpass(real(txDataUpsampled), 0.85*fs, fsCh, ImpulseResponse="iir", Steepness=0.8);
    txDataUpsampledImag = lowpass(imag(txDataUpsampled), 0.85*fs, fsCh, ImpulseResponse="iir", Steepness=0.8);
    txDataUpsampled     = txDataUpsampledReal + 1j*txDataUpsampledImag;
  end

  % Time vectors for upconversion and downconversion
  tDuration = length(txDataUpsampled) / fsCh;
  tup       = reshape(0 : 1/fsCh : (tDuration-1/fsCh), [], 1);

  mixerUp   = exp(1i*2*pi*fc*tup);
  mixerDown = exp(-1i*2*pi*fc*tup);

  % Upconversion mixing
  txDataMixUp = txDataUpsampled .* mixerUp;

  if strcmpi(modelType, 'rf')
    txDataMixUp = real(txDataMixUp); % When rf is used, we need to take the real part of the signal as we emulate passband signals
  end

  % Pass through the channel
  if strcmpi(channelType, 'rayleigh') || strcmpi(channelType, 'rician')
    txCh = channelObj(txDataMixUp);
  else
    txCh = txDataMixUp;
  end

  txDownMixed = txCh .* mixerDown;

  if strcmpi(modelType, 'rf')
    % Downconversion mixing
    txDownMixedReal = lowpass(real(txDownMixed), 0.85*fs, fsCh, ImpulseResponse="iir", Steepness=0.8);
    txDownMixedImag = lowpass(imag(txDownMixed), 0.85*fs, fsCh, ImpulseResponse="iir", Steepness=0.8);
    txDownMixed     = txDownMixedReal + 1j*txDownMixedImag;
  end

  txDownSampled = resample(txDownMixed, fs, fsCh);

  if randomPhaseOffset
    phaseOffset   = rand() * 2 * pi;
    txDownSampled = txDownSampled * exp(1i * phaseOffset);
  end

  if noiseEnable
    txDownSampled = awgn(txDownSampled, snrdB, 'measured');
  end

  rxDataArray(iFc, :) = txDownSampled;

  txHtltf    = txDownSampled(ind.HTLTF(1) : ind.HTLTF(2), :);
  txHtltfSym = wlanHTLTFDemodulate(txHtltf, cfgHT, symOffset);
  csiHTLTF   = wlanHTLTFChannelEstimate(txHtltfSym, cfgHT, smoothingSpan);

  csiSegmentsArray(iFc, :) = csiHTLTF;
end

% REVISIT: Update to match doppler in naming
[csiSpliced, csiSegmentsSpliced] = spliceCombineCSI(infoSplice, csiSegmentsArray);


%% HTLTF results
%------------------------------------------------------------------------------

if plotAllCFR
  nPlots = size(csiSegmentsArray,1);
else
  nPlots = 1;
end

for i = 1:nPlots
  csi = reshape(csiSegmentsArray(i, :), [], 1);
  cir = prunedIDFTMatHTLTF * csi;

  if strcmpi(modelType, 'rf')
    freqOffset   = infoSplice.CarrierFrequenciesPassband(i) / 1e9;
    freqMHzHTLTF = infoHTLTF.ActiveFrequencyIndices * scSpacing/1e9 + freqOffset;
  else
    freqOffset   = infoSplice.CarrierFrequenciesBaseband(i) / 1e6;
    freqMHzHTLTF = infoHTLTF.ActiveFrequencyIndices * scSpacing/1e6 + freqOffset;
  end
  timeSecHTLTF = (0:1:(length(cir)-1))/fs;

  plotComplexResponse(freqMHzHTLTF, csi, 'title', 'Est. CFR HT-LTF', 'xLabel', 'Freq. (MHz)', ...
    'useDbForMag', false, 'unwrapPhase', true, 'markerStyle', '-o', 'useStem', false);
  plotComplexResponse(timeSecHTLTF, cir, 'title', 'Est.s CIR HT-LTF', 'xLabel', 'Time (s)', ...
    'useDbForMag', false, 'unwrapPhase', false, 'markerStyle', '-o', 'useStem', true);
end


%% Plot segments
%-------------------------------------------------------------------------------
plotCSISegments(csiSegmentsArray, infoSplice, 'Unspliced segments'); % Show uncorrected segments
plotCSISegments(csiSegmentsSpliced, infoSplice, 'Spliced segments'); % Show corrected segments

%% Splicing results
%-------------------------------------------------------------------------------
freqMHzSplice  = infoSplice.ActiveFrequencyIndices * scSpacing/1e6;
timeSecSpliced = transpose((0:1:(infoSplice.FFTLength-1))/ infoSplice.SampleRate);

csiSpliced = csiSpliced .* exp(-1j*2*pi*freqMHzSplice*0.025); % TODO: Hack, find out what we put in here to shift the impulse response right

csiSplicedFilt = filterCSI(infoSplice, csiSpliced, smoothingWindowSize);
%cirSpliced     = ifftshift(prunedIDFTMatSpliced * transpose(csiSplicedFilt), 2); % ifftshift can help when phase is not exact and your main peak is shifting around
cirSpliced     = prunedIDFTMatSpliced * transpose(csiSplicedFilt);


plotComplexResponse(freqMHzSplice, csiSpliced, 'title', 'Est. CFR SPLICING', 'xLabel', 'Freq. (MHz)', ...
  'useDbForMag', true, 'unwrapPhase', true, 'markerStyle', '-o', 'useStem', false);
plotComplexResponse(timeSecSpliced, cirSpliced, 'title', 'Est.s CIR SPLICING', 'xLabel', 'Time (s)', ...
  'useDbForMag', false, 'unwrapPhase', false, 'markerStyle', '-o', 'useStem', true);

% Print the dimensions of csiSpliced
disp('Dimensions of csiSpliced:');
disp(size(csiSpliced));

% Print the dimensions of cirSpliced
disp('Dimensions of cirSpliced:');
disp(size(cirSpliced));

% Perform IFFT to reconstruct CIR from CSI
reconstructedCIR = ifft(csiSpliced, [], 1);

% Plot the original CIR
figure;
subplot(2, 1, 1);
stem(abs(cirSpliced), 'b', 'DisplayName', 'Original CIR', 'MarkerFaceColor', 'b');
title('Original CIR');
xlabel('Sample Index');
ylabel('Magnitude');
legend;
grid on;

% Plot the reconstructed CIR
subplot(2, 1, 2);
plot(abs(reconstructedCIR), 'r', 'DisplayName', 'Reconstructed CIR');
title('Reconstructed CIR from CSI');
xlabel('Sample Index');
ylabel('Magnitude');
legend;
grid on;


sgtitle('Comparison of Original and Reconstructed CIR');

set(gcf, 'Position', [100, 100, 800, 600]);

%% Save the data 
%-------------------------------------------------------------------------------

save('test_splicing.mat', '-v7.3');

%% True channel
%-------------------------------------------------------------------------------

% Conjugated version to match MATLAB visualizer
visualizeChannel(channelObj, fsSplicedChannel, true, plotMatlab);
if plotMatlab
  plotComplexResponse(-freqMHzSplice, conj(csiSpliced), 'title', 'Est. CFR SPLICING (X*(-jw))', 'xLabel', 'Freq. (MHz)', ...
  'useDbForMag', true, 'unwrapPhase', true, 'markerStyle', '-o', 'useStem', false);
end
