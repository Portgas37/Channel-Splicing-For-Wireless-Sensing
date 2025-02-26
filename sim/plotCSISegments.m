% plotCSISegments - This function plots the CSI segments that are to be spliced (they have not been spliced yet). Plotting
%   pre-splicing can be useful to determine the best splicing points and for debugging splicing functions
%
% USAGE
%   plotCSISegments(csiSegments, infoSplice)
%
% INPUT PARAMETERS
%   csiSegments: Array containing the CSI segments (either spliced or unspliced). Shape is [nFc, nSc] with nFc being the
%     number of frequency channels and nSc the number of CSI samples.
%   infoSplice:  Structure with splicing information. Similat to the HTL info structure
%
% DETAILS
%   This function plots the magnitude and phase of each CSI segment provided.
%
function plotCSISegments(csiSegments, infoSplice, plotTitle)

  if nargin < 3
    plotTitle = '';
  end

  figure;

  nFc      = size(csiSegments, 1);
  colormap = lines(nFc); % This will create a colormap with nFc distinct colors

  xStart = infoSplice.CSIMapping{1}(1);
  xEnd   = infoSplice.CSIMapping{end}(end);

  % Plotting the magnitude
  subplot(2, 1, 1);
  hold on;
  for iFc = 1:nFc
    csi    = csiSegments(iFc, :);
    csiMag = abs(csi);

    plot(infoSplice.CSIMapping{iFc}, csiMag, 'Color', colormap(iFc, :));

    % Draw a vertical line at the start of each segment and label it with the segment index
    xline(infoSplice.CSIMapping{iFc}(1),'-',sprintf('%d', iFc), 'Color', colormap(iFc, :), 'LabelOrientation', 'horizontal', 'LineWidth', 2);
  end
  xlabel('Spliced CSI index');
  ylabel('Magnitude');
  xlim([xStart, xEnd]);
  title([plotTitle ' (Magnitude)']);
  hold off;

  % Plotting the phase
  subplot(2, 1, 2);
  hold on;
  for iFc = 1:nFc
    csi      = csiSegments(iFc, :);
    csiPhase = unwrap(angle(csi));

    plot(infoSplice.CSIMapping{iFc}, csiPhase, 'Color', colormap(iFc, :));

    % Draw a vertical line at the start of each segment and label it with the segment index
    xline(infoSplice.CSIMapping{iFc}(1),'-',sprintf('%d', iFc), 'Color', colormap(iFc, :), 'LabelOrientation', 'horizontal', 'LineWidth', 2);
  end
  xlabel('Spliced CSI index');
  ylabel('Phase (radians)');
  xlim([xStart, xEnd]);
  title([plotTitle ' (Phase)']);
  hold off;

end
