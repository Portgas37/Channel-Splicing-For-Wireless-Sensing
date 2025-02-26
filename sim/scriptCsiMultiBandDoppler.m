% scriptCsiMultiBandDoppler.m.m
%
% This script aims to demonstrate the process of combining Channel State Information (CSI)
% from multiple carriers to perform channel splicing. It builds on the simpler 'scriptCsiSingle.m'
%
% PERFORMANCE
%   - Note that if you use rf modelType, the performance will be worse than with bb. This is because
%     we have to simulate at GHz frequencies.
%
%   - When you use complex channels, the performance will be much worse due to the Doppler simulation.
%
% REFERENCES
%   - Interpolation/Decimation: Figure 11.12 in Software Defined Radio using MATLAB & Simulink and the RTL-SDR

%% Settings
%------------------------------------------------------------------------------

clear; clc;
close all;
rng("default");

run('spliceSettingsCommon.m');

% REVISIT: Put all these into a cfg. so people can see it for the data

% Parameters
scenarioType      = 'real'; % One of ('test' | 'real')
scenarioIndex     = 1;
noiseEnable       = true;
channelType       = 'rayleigh'; % One of ('rayleigh' | 'rician' | 'none')
randomPhaseOffset = false; % Add a random phase offset to the signal. This represents the random phase offset
                          % That occurs when we change the carrier frequency.

snrdB               = 30;
symOffset           = 1.0;
smoothingSpan       = 1; % MATLAB smoothing function
smoothingWindowSize = 5; % Additional smoothing

modelType = 'bb'; % One of ('bb' | 'rf'). If 'bb', the signal stays in baseband, if rf, the signal is mixed up to a carrier frequency and converted to a passband signal

plotAllCFR = false; % ( false | true ). If true, plot all intermediate CFRs used in the splicing process
plotMatlab = false; % ( false | true ). If true, plot the channel as MATLAB visualizes it. This is useful to compare the results.
                    % However, please read the documentation of visualizeChannel() to understand the limitations of this function.

% Channel hopping parameters
fc0    = 2.4e9 + 20e6; % First carrier frequency
fcStep = 32;           % Jump X carrier spacings. If fcStep=16 jump 5 MHz which is Wi-Fi channel spacing
nFc    = 7;           % Number of channels to hop

fcList = fc0 + (0:(nFc-1)) * fcStep * scSpacing;

% Doppler parameters
fDoppler = 0.1; % Max Doppler frequency shift in Hz
nFrames  = 200;  % Number of frames to simulate


%% Splicing setup
%------------------------------------------------------------------------------

infoSplice = spliceInitInfo(infoHTLTF, scSpacing, fcList);
infoHTLTF.fftDcIdx   = getFftDcIdx(infoHTLTF.FFTLength);
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
  channelObj = generateChannel(channelType, scenario, fsCh, fDoppler);
end

%% Simulate
%------------------------------------------------------------------------------

txData = wlanWaveformGenerator([], cfgHT);

txDataArray = zeros(nFrames, nFc, length(txData), 'like', 1j);
rxDataArray = zeros(nFrames, nFc, length(txData), 'like', 1j);

csiSegmentsArray        = zeros(nFrames, nFc, infoSplice.infoWLANField.NumTones, 'like', 1j);
csiSegmentsSplicedArray = zeros(nFrames, nFc, infoSplice.infoWLANField.NumTones, 'like', 1j);
csiSplicedArray         = zeros(nFrames, infoSplice.NumTones, 'like', 1j);

txData = wlanWaveformGenerator([], cfgHT);
 
timeStepSize  = 20*length(txData)/fs
initTime      = 0;
timeMax       = timeStepSize*nFrames*nFc - timeStepSize;
timeStepArray = 0:timeStepSize:timeMax;

simStep = 0;

for iFrame = 1:nFrames
  for iFc = 1:nFc
    simStep = simStep + 1;
    txDataArray(iFrame, iFc, :) = txData;

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
      txCh = channelObj(txDataMixUp, timeStepArray(simStep)); % REVISIT: MOve out of here, so no timestamp in jump,, I just need different timesteps for inside a frame and outside
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

    rxDataArray(iFrame, iFc, :) = txDownSampled;

    txHtltf    = txDownSampled(ind.HTLTF(1) : ind.HTLTF(2), :);
    txHtltfSym = wlanHTLTFDemodulate(txHtltf, cfgHT, symOffset);
    csiHTLTF   = wlanHTLTFChannelEstimate(txHtltfSym, cfgHT, smoothingSpan);

    csiSegmentsArray(iFrame, iFc, :) = csiHTLTF;
  end

  [csiSpliced, csiSegmentsSpliced] = spliceCombineCSI(infoSplice, squeeze(csiSegmentsArray(iFrame, :, :)));

  csiSegmentsSplicedArray(iFrame, :, :) = csiSegmentsSpliced;
  csiSplicedArray(iFrame, :)            = csiSpliced;
end


%% HTLTF results
%------------------------------------------------------------------------------

if plotAllCFR
  nPlots = size(csiSegmentsArray,2);
else
  nPlots = 1;
end

frameToPlot = 1;

for i = 1:nPlots
  csi = reshape(csiSegmentsArray(frameToPlot, i, :), [], 1);
  cir = prunedIDFTMatHTLTF * csi;

  if strcmpi(modelType, 'rf')
    plotFreqOffset = infoSplice.CarrierFrequenciesPassband(i) / 1e9;
    freqMHzHTLTF   = infoHTLTF.ActiveFrequencyIndices * scSpacing/1e9 + plotFreqOffset;
  else
    plotFreqOffset = infoSplice.CarrierFrequenciesBaseband(i) / 1e6;
    freqMHzHTLTF   = infoHTLTF.ActiveFrequencyIndices * scSpacing/1e6 + plotFreqOffset;
  end
  timeSecHTLTF = (0:1:(length(cir)-1))/fs;

  plotComplexResponse(freqMHzHTLTF, csi, 'title', 'Est. CFR HT-LTF', 'xLabel', 'Freq. (MHz)', ...
    'useDbForMag', false, 'unwrapPhase', true, 'markerStyle', '-o', 'useStem', false);
  plotComplexResponse(timeSecHTLTF, cir, 'title', 'Est.s CIR HT-LTF', 'xLabel', 'Time (s)', ...
    'useDbForMag', false, 'unwrapPhase', false, 'markerStyle', '-o', 'useStem', true);
end


%% Plot segments
%-------------------------------------------------------------------------------

% Just show the first and last frames
plotCSISegments(squeeze(csiSegmentsArray(1, :, :)), infoSplice, 'Unspliced segments (first frame)');      % Show uncorrected segments
plotCSISegments(squeeze(csiSegmentsSplicedArray(1, :, :)), infoSplice, 'Spliced segments (first frame)'); % Show corrected segments

plotCSISegments(squeeze(csiSegmentsArray(end, :, :)), infoSplice, 'Unspliced segments (last frame)');
plotCSISegments(squeeze(csiSegmentsSplicedArray(end, :, :)), infoSplice, 'Spliced segments (last frame)');

%% Splicing results
%-------------------------------------------------------------------------------

plotIndices    = [1, size(csiSplicedArray,1)]; % These are all the frames you want to plot
freqMHzSplice  = infoSplice.ActiveFrequencyIndices * scSpacing/1e6;
timeSecSpliced = transpose((0:1:(infoSplice.FFTLength-1))/ infoSplice.SampleRate);

for i = 1:length(plotIndices)
  plotIdx = plotIndices(i);
  plotStr = sprintf('(Tstep=%.3f sec)', timeStepArray((plotIdx-1)*nFc+1));

  csiSplicedFilt = filterCSI(infoSplice, csiSplicedArray(plotIdx, :), smoothingWindowSize);
  cirSpliced     = prunedIDFTMatSpliced * transpose(csiSplicedFilt);

  titleCFR = sprintf('Est. CFR SPLICING %s', plotStr);
  titleCIR = sprintf('Est. CIR SPLICING %s', plotStr);

  plotComplexResponse(freqMHzSplice, csiSplicedArray(plotIdx, :), 'title', titleCFR, ...
    'xLabel', 'Freq. (MHz)', 'useDbForMag', true, 'unwrapPhase', true, 'markerStyle', '-o', 'useStem', false);
  plotComplexResponse(timeSecSpliced, cirSpliced, 'title', titleCIR, ...
    'xLabel', 'Time (s)', 'useDbForMag', false, 'unwrapPhase', false, 'markerStyle', '-o', 'useStem', true);
end

%% Save the data 
%-------------------------------------------------------------------------------

save('test_splicing_doppler.mat', '-v7.3');


%% True channel
%-------------------------------------------------------------------------------

% Conjugated version to match MATLAB visualizer
visualizeChannel(channelObj, fsSplicedChannel, true, plotMatlab, timeStepArray(end-nFc+1));
if plotMatlab
  plotComplexResponse(-freqMHzSplice, conj(csiSpliced), 'title', 'Est. CFR SPLICING (X*(-jw))', 'xLabel', 'Freq. (MHz)', ...
    'useDbForMag', true, 'unwrapPhase', true, 'markerStyle', '-o', 'useStem', false);
end
