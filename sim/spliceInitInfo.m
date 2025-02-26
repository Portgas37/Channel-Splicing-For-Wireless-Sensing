%spliceInitInfo: This function initializes the splicing information.
%
% USAGE
%   infoSplice = spliceInitInfo(infoWLANField, scSpacing, fcList)
%
% INPUT PARAMETERS
%   infoWLANField : WLAN field information. Either for  HT Long Training Field (HT-LTF), nonHT etc.
%   scSpacing     : Subcarrier spacing.
%   fcList        : List of carrier frequencies (passband)
%
% OUTPUT PARAMETERS
%   infoSplice                    : A structure containing the initialized splicing information.
%     .ActiveFrequencies          : Frequencies in the baseband.
%     .ActiveFrequenciesPassband  : Frequencies in the passband.
%     .ActiveFFTIndices           : Indices into the FFT matrix.
%     .BwIncrease                 : Increase in bandwidth due to splicing.
%     .CarrierFrequenciesBaseband : Baseband frequencies for each carrier.
%     .CarrierFrequenciesPassBand : Passband frequencies for each carrier.
%     .CSIMapping                 : Mapping into the new structure.
%     .FFTLength                  : FFT length after splicing.
%     .FFTIncrease                : Increment in FFT length due to splicing.
%     .fftDcIdx                   : Index position of the DC component in the FFT.
%     .infoWLANField              : WLAN field information from which we spliced from at each carrier
%     .NumTones                   : Number of tones (or subcarriers) after splicing.
%     .SampleRate                 : Sample rate after splicing.
%
% DETAILS
%   This function provides an initialization for the splicing information, which is equivalent to
%   the wlanHTOFDMInfo structure but adapted for spliced carriers.

function infoSplice = spliceInitInfo(infoWLANField, scSpacing, fcList)
  fcList = sort(fcList);
  fc0    = fcList(1) + (fcList(end) - fcList(1)) / 2; % fc0 is the frequency in the middle based on the end points

  % Check that fc0 is a multiple of scSpacing
  if mod(fc0, scSpacing) > 1e-10
    error('fc0 is not a multiple of scSpacing!');
  end

  infoSplice = initializeSpliceStruct(infoWLANField);
  infoSplice = updateSpliceInfo(infoSplice, infoWLANField, fcList, scSpacing, fc0);
  infoSplice.ActiveFrequenciesPassband = infoSplice.ActiveFrequencies + fc0;
  infoSplice.fc0 = fc0;

  if ~isequal(infoSplice.CarrierFrequenciesPassband(:), fcList(:))
    warning('infoSplice.CarrierFrequenciesPassband and fcList are not identical.');
  end

  fprintf("Splicing information:\n");
  fprintf("\tOriginal bandwidth = %d MHz\n", infoWLANField.SampleRate / 1e6);
  fprintf("\tNew bandwidth = %d MHz\n", infoSplice.SampleRate / 1e6);
end


function infoSplice = initializeSpliceStruct(infoWLANField)
  infoSplice = struct();
  infoSplice.ActiveFrequencies          = [];
  infoSplice.ActiveFrequenciesPassband  = [];
  infoSplice.ActiveFFTIndices           = [];
  infoSplice.BwIncrease                 = NaN;
  infoSplice.CarrierFrequenciesBaseband = [];
  infoSplice.CSIMapping                 = {};
  infoSplice.FFTIncrease                = NaN;
  infoSplice.FFTLength                  = NaN;
  infoSplice.NumTones                   = NaN;
  infoSplice.SampleRate                 = NaN;
  infoSplice.fftDcIdx                   = NaN;

  infoSplice.infoWLANField = infoWLANField; % WLAN field information from which we spliced from
end


function infoSplice = updateSpliceInfo(infoSplice, infoWLANField, fcList, scSpacing, fc0)

  % Each channel has the same bandwidth and FFT size, so the FFT increase is based on assuming a full channel from left to right based on fcList
  spliceFFTSize = (fcList(end)-fcList(1))/scSpacing + infoWLANField.FFTLength;

  infoSplice.FFTIncrease = spliceFFTSize - infoWLANField.FFTLength;
  infoSplice.BwIncrease  = infoSplice.FFTIncrease * scSpacing;
  infoSplice.FFTLength   = spliceFFTSize;
  infoSplice.SampleRate  = infoWLANField.SampleRate + infoSplice.BwIncrease;
  infoSplice.fftDcIdx    = getFftDcIdx(infoSplice.FFTLength);

  % REVISIT: Ensure shapes
  infoSplice = computeFFTIndicesForSplicing(infoSplice, infoWLANField, fcList, scSpacing);
  infoSplice.ActiveFrequencyIndices     = infoSplice.ActiveFFTIndices - infoSplice.fftDcIdx;
  infoSplice.CarrierFrequenciesPassband = infoSplice.CarrierFrequenciesBaseband + fc0;

  infoSplice.NumTones = length(infoSplice.ActiveFrequencies);
end


%computeFFTIndicesForSplicing - Computes the FFT indices for spliced carriers and determines their active frequencies.
%
% USAGE
%   infoSplice = computeFFTIndicesForSplicing(infoSplice, infoWLANField, nFc, fcStep)
%
% INPUT PARAMETERS
%   infoSplice    - Current splicing information.
%   infoWLANField - Information about WLAN field that we base the splicing on.
%   nFc           - Number of spliced channels.
%   fcStep        - Step in terms of sub-carriers between carriers.
%   scSpacing     - Subcarrier spacing.
%
% OUTPUT PARAMETERS
%   infoSplice - Updated splicing information with new FFT indices and other related data.
%
% DETAILS
%   Calculates FFT indices based on carrier positions relative to reference frequency (fc0), determining
%   appropriate active FFT indices and the DC component's position. The function also returns a CSI mapping
%   for each carrier.
%
function infoSplice = computeFFTIndicesForSplicing(infoSplice, infoWLANField, fcList, scSpacing)
  freqVec = generateFreqAxis(infoSplice.SampleRate, infoSplice.FFTLength);

  for iFc = 1:length(fcList)
    fftIdxOffset = (fcList(iFc) - fcList(1))/scSpacing;

    activeFFTIndices            = fftIdxOffset + infoWLANField.ActiveFFTIndices;
    infoSplice.ActiveFFTIndices = union(infoSplice.ActiveFFTIndices, activeFFTIndices);
    infoSplice.CSIMapping{iFc}  = activeFFTIndices;

    carrierOffset       = fftIdxOffset + getFftDcIdx(infoWLANField.FFTLength);
    carrierFreqBaseband = freqVec(carrierOffset);
    infoSplice.CarrierFrequenciesBaseband(iFc) = carrierFreqBaseband;

    activeFreq = carrierFreqBaseband + infoWLANField.ActiveFrequencyIndices * scSpacing;
    infoSplice.ActiveFrequencies = union(infoSplice.ActiveFrequencies, activeFreq);
  end
end

function fftDcIdx = getFftDcIdx(nfft)
  fftDcIdx = floor(nfft / 2) + 1;
end
