% generateChannel: This function generates a wireless channel model, either Rayleigh or Rician.
%
% USAGE
%   channelObj = generateChannel(channelType, selectedScenario, fsChannel)
%
% INPUT PARAMETERS
%   channelType       : Type of the channel ('rayleigh' or 'rician')
%   selectedScenario  : Scenario struct containing 'delays' and 'attenuations' for the channel
%   fsChannel         : Sample rate of the channel [Hz] after mixing to carrier (just the signal BW if in baseband)
%   fDoppler          : Maximum Doppler shift [Hz] (default: 0.0). For a person walking, it would be around 10-20 Hz
%
% OUTPUT PARAMETERS
%   channelObj: An object of the generated channel
%
% DETAILS
%   The user can define specific scenarios with path delays and average path gains.
%   When you mix to fc you get that the bandwidth of the channel is fsChannel.
%
% REFERENCE
% - <https://www.mathworks.com/help/comm/ref/comm.rayleighchannel-system-object.html>
%
function channelObj = generateChannel(channelType, selectedScenario, fsChannel, fDoppler)

  if nargin < 4
    fDoppler = 0.0;
  end

  % Note that when 'Visualization' is enabled, the channel is plotted as soon as you apply the channel
  % to some data). If the Doppler is 0, we get a warning when using the InitialTimeSource as 'Input port'.
  % Therefore, we only enable the visualization in very specific cases like in the visualizeChannel() function.
  settingsRayleigh = {
    'SampleRate', fsChannel, ...
    'PathDelays', selectedScenario.delays, ...
    'AveragePathGains', selectedScenario.attenuations, ...
    'NormalizePathGains', true, ...
    'MaximumDopplerShift', fDoppler, ...
    'RandomStream', 'mt19937ar with seed', ...
    'Seed', 0, ...
  };

  settingsRician = {
    'KFactor', 3.0, ...
    'SampleRate', fsChannel, ...
    'PathDelays', selectedScenario.delays, ...
    'AveragePathGains', selectedScenario.attenuations, ...
    'NormalizePathGains', true, ...
    'MaximumDopplerShift', fDoppler, ...
    'RandomStream', 'mt19937ar with seed', ...
    'Seed', 0, ...
  };

  % Create the channel based on type
  switch lower(channelType)
    case 'rayleigh'
      channelObj = comm.RayleighChannel(settingsRayleigh{:});
    case 'rician'
      channelObj = comm.RicianChannel(settingsRician{:});
    otherwise
      error('Unsupported channel type.');
  end

  if fDoppler > 0
    channelObj.InitialTimeSource = 'Input port';
    channelObj.FadingTechnique   = 'Sum of sinusoids';
  end
end
