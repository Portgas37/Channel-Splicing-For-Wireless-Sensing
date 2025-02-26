% csiSplicerPhases - Splices phases from multiple CSI estimates together.
%
% USAGE
%   [splicedMag, splicedPhase] = csiSplicerPhases(infoSplice, csi2Splice, overlapSize)
%
% INPUT PARAMETERS
%   infoSplice - Current splicing information.
%   csi2Splice - Matrix containing CSI data from multiple carriers.
%   overlapSize - Size of boundary fix region. Default is 3.
%
% OUTPUT PARAMETERS
%   splicedMag   - Magnitude information after splicing.
%   splicedPhase - Phase information after splicing.
%
% DESCRIPTION
%   This function splices together CSI data from multiple carriers by aligning their phases.
%   Since segments have to overlap for this method to be spliced together, we use the fact that
%   the phase is continuous and use the median phase difference between overlapping segments to
%   correct for the random phase offset between segments.
%
%   We also correct the boundaries of each segment using values from the previous or next segment.
%
%   Note that a phase slope cannot be used for the offset due to how the phase oscillates in a
%   multi-path environment. I have some old code at the end of this file that uses a phase slope
%   just to show that it doesn't work.

function [splicedMag, splicedPhase] = csiSplicerPhases(infoSplice, csi2Splice, endPointSize)

  if nargin < 3 || isempty(endPointSize)
    endPointSize = 3;
  end

  csi2SpliceSize = size(csi2Splice);
  nFc = csi2SpliceSize(1);
  nSc = csi2SpliceSize(end);

  % For unwrapping, just calculate difference
  splicedPhase = unwrap(angle(csi2Splice), [], 2);
  splicedMag   = abs(csi2Splice);

  for idxFc = 1:(nFc-1)
    currMapping = infoSplice.CSIMapping{idxFc};
    nextMapping = infoSplice.CSIMapping{idxFc+1};

    % Find pointers to where they overlap and use these for splicing
    overlapMapping = intersect(currMapping, nextMapping);
    overlapSize    = length(overlapMapping);

    % Find the indices in currMapping that correspond to overlapMapping, i.e., currMapping(overlapIndices) = overlapMapping
    [~, currIndices] = ismember(overlapMapping, currMapping);
    % Find the indices in nextMapping that correspond to overlapMapping, i.e., nextMapping(overlapIndices) = overlapMapping
    [~, nextIndices] = ismember(overlapMapping, nextMapping);

    % Compute phase difference between current and next segment
    overlapPhaseCurr = splicedPhase(idxFc, currIndices);
    overlapPhaseNext = splicedPhase(idxFc+1, nextIndices);

    phaseDiff = median(overlapPhaseCurr - overlapPhaseNext); % Median is robust to errors at edges etc.

    % Adjust the phase of the next segment based on the phase difference
    splicedPhase(idxFc+1, :) = splicedPhase(idxFc+1, :) + phaseDiff;

    % Minor correction at boundary, we correct boundaries of both
    startPointsCurr = currIndices(1:endPointSize);
    startPointsNext = nextIndices(1:endPointSize);
    endPointsCurr = currIndices(end - endPointSize + 1 : end);
    endPointsNext = nextIndices(end - endPointSize + 1 : end);

    % Correct current endpoints
    splicedPhase(idxFc, endPointsCurr) = splicedPhase(idxFc+1, endPointsNext);
    splicedMag(idxFc, endPointsCurr)   = splicedMag(idxFc+1, endPointsNext);

    % Correct next start points
    splicedPhase(idxFc+1, startPointsNext) = splicedPhase(idxFc, startPointsCurr);
    splicedMag(idxFc+1, startPointsNext)   = splicedMag(idxFc, startPointsCurr);
  end

  splicedPhase = splicedPhase - splicedPhase(1); % Normalize phase so phase is 0 at first subcarrier
end

% iFc = 5; figure(1); plot(infoSplice.CSIMapping{iFc}, splicedPhase(iFc, :)); hold on; plot(infoSplice.CSIMapping{iFc+1}, splicedPhase(iFc+1, :)); hold off;
% iFc = 5; figure(2); plot(infoSplice.CSIMapping{iFc}, splicedMag(iFc, :)); hold on; plot(infoSplice.CSIMapping{iFc+1}, splicedMag(iFc+1, :)); hold off;
% iFc = 5; figure(3); plot(infoSplice.CSIMapping{iFc}, 20*log10(splicedMag(iFc, :))); hold on; plot(infoSplice.CSIMapping{iFc+1}, 20*log10(splicedMag(iFc+1, :))); hold off;








%

%
% function [splicedMag, splicedPhase] = csiSplicerPhases(infoSplice, csi2Splice, overlapSize)
%   if nargin < 3 || isempty(overlapSize)
%     overlapSize = 3;
%   end
%
%   csi2SpliceSize = size(csi2Splice);
%   nFc = csi2SpliceSize(1);
%   nSc = csi2SpliceSize(end);
%
%   % For unwrapping, just calculate difference
%   splicedPhase = unwrap(angle(csi2Splice), [], 2);
%   splicedMag   = abs(csi2Splice);
%
%   % Fit a line to the phase of the first CSI for reference
%   firstCsiMapping = infoSplice.CSIMapping{1};
%   firstCsiPhase   = splicedPhase(1, :);  % We don't index with firstCsiMapping here
%   P1 = polyfit(firstCsiMapping, firstCsiPhase, 1); % First order polynomial fit
%   refLine = polyval(P1, firstCsiMapping);
%
%   for iFc = 2:nFc
%     % Fit a line to the phase of current CSI's inner subcarriers
%     currMapping  = infoSplice.CSIMapping{iFc};
%     currCsiPhase = splicedPhase(iFc, :);
%     P2 = polyfit(currMapping, currCsiPhase, 1);
%     currLine = polyval(P2, currMapping);
%
%     % Compute the phase offset between the reference line and the current line
%     % This offset provides a correction value ensuring the current segment aligns with the reference.
%     phaseOffset = mean(currLine - refLine);
%
%     % Adjust the phase of the current segment based on the phase offset
%     splicedPhase(iFc, :) = splicedPhase(iFc, :) - phaseOffset;
%   end
%
%   % Splicing loop for phase alignment
%   for iFc = 1:(nFc-1)
%     currMapping = infoSplice.CSIMapping{iFc};
%     nextMapping = infoSplice.CSIMapping{iFc+1};
%
%     % For fixing backward, find where the current mapping ends in the next mapping
%     nextIdxStart = nextMapping(end) - currMapping(end);
%
%     % Check if the mapping is a continuous block for both current and next mapping
%     diffCurrMapping = diff(currMapping);
%     diffNextMapping = diff(nextMapping);
%
%     if length(unique(diffCurrMapping)) == 1 && length(unique(diffNextMapping)) == 1
%       splicedPhase(iFc, end - overlapSize + 1:end) = splicedPhase(iFc+1, nextIdxStart : nextIdxStart + overlapSize - 1);
%       splicedMag(iFc, end - overlapSize + 1:end)   = splicedMag(iFc+1, nextIdxStart : nextIdxStart + overlapSize - 1);
%     end
%   end
%
%   % Normalize the phase such that the phase is 0 at the first subcarrier
%   splicedPhase = splicedPhase - splicedPhase(1);
% end
