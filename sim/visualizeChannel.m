% visualizeChannel: This function visualizes the generated wireless channel model from MATLAB.
%
% USAGE
%   visualizeChannel(channelObj, fsChannel)
%   visualizeChannel(channelObj, fsChannel, plotEn)
%   visualizeChannel(channelObj, fsChannel, plotEn, matlabPlotEn)
%
% INPUT PARAMETERS
%   channelObj   : An object of the generated channel
%   fsChannel    : Sample rate (bandwidth) of the channel [Hz] we want to measure with splicing
%   plotEn       : Boolean, true if channel plots (impulse and frequency response) are desired
%   matlabPlotEn : Boolean, true if MATLAB's built-in visualization is enabled
%
% DETAILS
%   This function provides visualization for the impulse response and frequency response of the channel.
%
%   Note that the matlab plot is always based on the baseband signal, you might see that the response is
%   shifted if you try to simulate at passband.
%
% VISUALIZATION NOTES
%  - For some reason, when MATLAB visualizes the channel, it does match what we get from
%    sending an impulse response through the channel. To match the MATLAB visualizer we need to apply the relationship
%    x*(t) <--> X*(-jw) which can be achieved by plotting with negative frequency axis, i.e., -freqMHzSplice, instead
%    and conjugating the frequency response. In the time-domain, this will just change the phase of the taps.
%    but not the magnitude.
%
% - The MATLAB plotting does not work with Doppler turned on, only our own custom plotting works with Doppler.
%
function visualizeChannel(channelObj, fsChannel, plotEn, matlabPlotEn, currentTime)

  if nargin < 3 || isempty(plotEn), plotEn = false; end
  if nargin < 4 || isempty(matlabPlotEn), matlabPlotEn = false; end
  if nargin < 5, currentTime = 0; end

  if plotEn || matlabPlotEn
    release(channelObj); % Critical for correct MATLAB visualization and impulse response
    channelObj.SampleRate = fsChannel;

    scSpacing     = 312.5e3;
    nfftChannel   = fsChannel / scSpacing;
    impulseSignal = [1; zeros(nfftChannel-1,1)]; % Impulse signal

    % If Doppler shift is greater than 0, simulate over time and show result from first and last step
    if channelObj.MaximumDopplerShift > 0.0
      timeVector = [0.0, currentTime];
      for i = 1:length(timeVector)
        impulseResponse = simulateChannel(channelObj, impulseSignal, timeVector(i));
        plotResponse(fsChannel, nfftChannel, impulseResponse, plotEn, matlabPlotEn, [' at time ' num2str(timeVector(i)) 's']);
      end

    % No Doppler shift, single simulation
    else
      impulseResponse = channelObj(impulseSignal);
      plotResponse(fsChannel, nfftChannel, impulseResponse, plotEn, matlabPlotEn);
    end

    release(channelObj); % Release after plotting
    if matlabPlotEn
      setupMatlabVisualization(channelObj);
    end
  end
end

function impulseResponse = simulateChannel(channelObj, impulseSignal, currentTime)
  if nargin < 3
    impulseResponse = channelObj(impulseSignal);
  else
    reset(channelObj); % Reset the channel to currentTime
    impulseResponse = channelObj(impulseSignal, currentTime);
  end
end

function plotResponse(fsChannel, nfftChannel, impulseResponse, plotEn, matlabPlotEn, titleSuffix)
  if nargin < 6, titleSuffix = ''; end

  % Sub-function for plotting
  if plotEn
    freqResponse = fftshift(fft(impulseResponse));
    freqVecMHz   = generateFreqAxis(fsChannel, nfftChannel) / 1e6;
    timeVec      = (0 : 1 : nfftChannel-1) / fsChannel;

    % Frequency and time-domain plots
    plotComplexResponse(freqVecMHz, freqResponse, 'title', ['TRUE CFR' titleSuffix], 'xLabel', 'Frequency (MHz)', ...
      'useDbForMag', true, 'unwrapPhase', true, 'markerStyle', '-o', 'useStem', false);

    if matlabPlotEn
      % Plot with negative frequency axis for MATLAB visualization
      plotComplexResponse(-freqVecMHz, conj(freqResponse), 'title', ['TRUE CFR (X*(-jw))' titleSuffix], 'xLabel', 'Frequency (MHz)', ...
        'useDbForMag', true, 'unwrapPhase', true, 'markerStyle', '-o', 'useStem', false);
    end

    plotComplexResponse(timeVec, impulseResponse, 'title', ['TRUE CIR' titleSuffix], 'xLabel', 'Time (s)', ...
      'useDbForMag', false, 'unwrapPhase', false, 'markerStyle', '-o', 'useStem', true);
  end
end

function setupMatlabVisualization(channelObj)
  % Setup MATLAB built-in visualization properties
  channelObj.Visualization    = 'Impulse and frequency responses';
  channelObj.SamplesToDisplay = '100%';
  channelObj.ChannelFiltering = 0;

  % Invoke the visualization and then release the object
  channelObj();
  release(channelObj);
  channelObj.Visualization    = 'Off';
  channelObj.ChannelFiltering = 1;
end
