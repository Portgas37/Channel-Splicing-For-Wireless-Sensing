% filterCSI - Filters the spliced CSI data with a median filter, handling gaps. Note that large gaps
%   are perfectly fine as we just filter in the region where we have sub-carriers.
%
% USAGE
%   filtCSI = filterCSI(infoSplice, csi, windowSize)
%
% INPUT PARAMETERS
%   infoSplice  - Splicing information containing active frequency indices.
%   csi         - Spliced CSI information containing both magnitude and phase. Shape must be [nFrames, nSubcarriers].
%   windowSize  - (Optional) Maximum window size for median filtering, default is 3.
%
% OUTPUT PARAMETERS
%   filtCSI - Filtered CSI information
function [filtCSI, filtCSIFull] = filterCSI(infoSplice, csi, windowSize)
  if nargin < 3
    windowSize = 3;  % Default maximum window size for median filtering
  end

  if isvector(csi)
    csi = reshape(csi, 1, []);
  end

  activeIndices = infoSplice.ActiveFFTIndices;
  nFrames       = size(csi, 1);

  magFull   = NaN(nFrames, infoSplice.FFTLength);
  phaseFull = NaN(nFrames, infoSplice.FFTLength);
  magFull(:, activeIndices)   = abs(csi);
  phaseFull(:, activeIndices) = unwrap(angle(csi), [], 2);

  % Fill end points and remaining with nearest
  magFull   = fillmissing(magFull, 'nearest', 2, 'EndValues', 'nearest');
  phaseFull = fillmissing(phaseFull, 'nearest', 2, 'EndValues', 'nearest');

  magFull   = movmean(magFull, windowSize, 2);
  phaseFull = movmean(phaseFull, windowSize, 2);
  %magFull   = movmedian(magFull, windowSize, 2);
  %phaseFull = movmedian(phaseFull, windowSize, 2);


  % Post-processing - combine magnitude and phase
  filtCSIFull = magFull .* exp(1j * phaseFull);
  filtCSI     = filtCSIFull(:, activeIndices);
end
