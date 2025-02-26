%plotComplexResponse: Plots magnitude and phase of complex-valued signals. This function plots the
%   magnitude and phase of complex-valued signals given the x-axis values and a matrix of complex
%   signals, where each column represents a separate signal and the rows are time.
%
%   USAGE:
%   plotComplexResponse(x, y)
%   plotComplexResponse(x, y, 'Parameter', value, ...)
%
%   INPUT PARAMETERS:
%   x: Vector of x-axis values.
%   y: 2D array of complex y-axis values (or vector), where each row is a signal.
%
%   OPTIONAL PARAMETERS:
%   'title'           : Main title for the plots (default 'Complex Signal Response').
%   'xLabel'          : Label for the x-axis (default 'x-values').
%   'yLabelMagnitude' : Label for the y-axis on the magnitude plot (default 'Magnitude').
%   'yLabelPhase'     : Label for the y-axis on the phase plot (default 'Phase (radians)').
%   'useDbForMag'     : Boolean, true to plot magnitude in dB (default true).
%   'unwrapPhase'     : Boolean, true to unwrap the phase (default false).
%   'markerStyle'     : Marker style for the plots (default '-o').
%   'useStem'         : Boolean, true to use stem plot instead of regular plot (default false).
%
%   EXAMPLES:
%   plotComplexResponse(frequency, response);
%   plotComplexResponse(time, signal, 'title', 'Signal Response', 'xLabel', 'Time (s)', 'useDbForMag', false);
%
function plotComplexResponse(x, y, varargin)
  % Define the parameters with defaults and add them to the input parser
  p = inputParser;
  addParameter(p, 'title', 'Complex Signal Response');
  addParameter(p, 'xLabel', 'x-values');
  addParameter(p, 'yLabelMagnitude', 'Magnitude');
  addParameter(p, 'yLabelPhase', 'Phase (radians)');
  addParameter(p, 'useDbForMag', true, @islogical);
  addParameter(p, 'unwrapPhase', false, @islogical);
  addParameter(p, 'markerStyle', '-o', @ischar);
  addParameter(p, 'useStem', false, @islogical);

  % Parse the input parameters
  parse(p, varargin{:});
  opts = p.Results;

  bIsVector = isvector(y);
  x = reshape(x, [], 1);

  % Pre-calculate the magnitude and phase if necessary
  if opts.useDbForMag
    magY     = 20 * log10(abs(y));
    magLabel = [opts.yLabelMagnitude ' (dB)'];
  else
    magY     = abs(y);
    magLabel = opts.yLabelMagnitude;
  end

  if opts.unwrapPhase
    if bIsVector
      phaseY = unwrap(angle(y));
    else
      phaseY = unwrap(angle(y), [], 2); % Unwrap along each row for a matrix (over time)
    end
  else
    phaseY = angle(y);
  end  

  % Create a new figure
  figure;

  % Plot magnitude
  subplot(2,1,1);
  plotSignalComponent(x, magY, opts, magLabel);
  title([opts.title ' (Magnitude)']);
  xlabel(opts.xLabel);
  grid on;

  % Plot phase
  subplot(2,1,2);
  plotSignalComponent(x, phaseY, opts, opts.yLabelPhase);
  title([opts.title ' (Phase)']);
  xlabel(opts.xLabel);
  grid on;
end


function plotSignalComponent(x, yComponent, opts, yLabel)
  hold on; % Hold on to plot multiple lines
  if isvector(yComponent)
    if opts.useStem
      stem(x, yComponent, opts.markerStyle);
    else
      plot(x, yComponent, opts.markerStyle);
    end
  else
    for k = 1:size(yComponent, 1) % Loop over each row of yComponent
      if opts.useStem
        stem(x, yComponent(k, :), opts.markerStyle);
      else
        plot(x, yComponent(k, :), opts.markerStyle);
      end
    end
  end
  hold off; % Release the hold after plotting all lines
  ylabel(yLabel);

  yMin = min(yComponent(:));
  yMax = max(yComponent(:));

  if abs(yMax - yMin) > 1e-6
    ylim([min(yComponent(:)), max(yComponent(:))]);
  end
end
