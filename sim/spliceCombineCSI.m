% spliceCombineCSI - Combines multiple CSI estimates by splicing them together to get a wider channel
%
% USAGE
%   combinedCSI = spliceCombineCSI(infoSplice, csi2Splice, method)
%
% INPUT PARAMETERS
%   infoSplice - Current splicing information.
%   csi2Splice   - Cell array containing CSI data from multiple carriers.
%   method      - String specifying the splicing method. Currently supported: 'splicePhases'.
%
% OUTPUT PARAMETERS
%   combinedCSI - Combined CSI information after splicing.
%
% DETAILS
%   The function combines CSI estimates from different carriers to create a spliced, composite CSI.
%   First, a splicing method is applied and afterwards the processed segments  are added together,
%   and the result is normalized.
%
%   The combination is driven by a splicing map provided in the infoSplice input, which dictates how
%   different CSI data are merged.
%
function [csiSpliced, csiSegments] = spliceCombineCSI(infoSplice, csi2Splice, method)

  if nargin < 3 || isempty(method)
    method = 'splicePhases';
  end

  csi2SpliceSize = size(csi2Splice);
  nFc = csi2SpliceSize(1);
  nSc = csi2SpliceSize(end);

  csiSplicedMag   = zeros(infoSplice.FFTLength, 1);
  csiSplicedPhase = zeros(infoSplice.FFTLength, 1);
  csiScaling      = zeros(infoSplice.FFTLength, 1, 'like', 1j);

  % Process the segments so they can be spliced together
  switch method
    case 'splicePhases'
      [csi2SpliceMag, csi2SplicePhase] = csiSplicerPhases(infoSplice, csi2Splice);
    otherwise
      error('Unsupported splicing method.');
  end

  % Post-processing: For methods that don't return the fully spliced data but need some common
  % post-processing. If your method returns the fully spliced data, you can just return at that point in the switch statement
  csiSegments = csi2SpliceMag .* exp(1j * csi2SplicePhase);

  % Combine the CSI estimates by adding them together (later divide so we average)
  for iFc = 1:nFc
    mapping = infoSplice.CSIMapping{iFc};
    csiSplicedMag(mapping)   = csiSplicedMag(mapping)   + transpose(csi2SpliceMag(iFc, :));
    csiSplicedPhase(mapping) = csiSplicedPhase(mapping) + transpose(csi2SplicePhase(iFc, :));
    csiScaling(mapping)      = csiScaling(mapping) + 1;
  end

  % Normalize the combined CSI
  csiSplicedMag   = csiSplicedMag ./ csiScaling;
  csiSplicedPhase = csiSplicedPhase ./ csiScaling;

  csiSpliced = csiSplicedMag .* exp(1j * csiSplicedPhase);
  csiSpliced = csiSpliced(infoSplice.ActiveFFTIndices); % Only take out the subset
end

% iFc = 5; figure(1); plot(infoSplice.CSIMapping{iFc}, unwrap(angle(csiSegments(iFc, :)))); hold on; plot(infoSplice.CSIMapping{iFc+1}, unwrap(angle(csiSegments(iFc+1, :)))); hold off
% iFc = 5; figure(2); plot(infoSplice.CSIMapping{iFc}, abs(csiSegments(iFc, :))); hold on; plot(infoSplice.CSIMapping{iFc+1}, abs(csiSegments(iFc+1, :))); hold off;
% iFc = 5; figure(3); plot(infoSplice.CSIMapping{iFc}, 20*log10(abs(csiSegments(iFc, :)))); hold on; plot(infoSplice.CSIMapping{iFc+1}, 20*log10(csiSegments(iFc+1, :))); hold off;
