%getMultipathScenario: Returns a multipath scenario given an index and computes the bandwidths required for each path.
%
% USAGE
%   [scenario, requiredBandwidth] = getMultipathScenario(scenarioType, scenarioIndex);
%   [scenario, requiredBandwidth] = getMultipathScenario(scenarioType, scenarioIndex, fs);
%
% INPUT PARAMETERS
%   scenarioType  : Type of the desired multipath scenario. Can be 'test' or 'real'.
%   scenarioIndex : Index of the desired multipath scenario.
%   fs            : Sampling frequency. Normally 20 MHz which is the bandwidth of 802.11a/g.
%
% OUTPUT PARAMETERS
%   scenario          : The selected multipath scenario structure.
%   requiredBandwidth : The bandwidths required for each path.
%
% DETAILS
%   This function provides different channel scenarios to illustrate various multipath effects.
%   It computes the individual bandwidths required to resolve each path, taking into consideration
%   the possible existence of a zero-delay path.
%
% NOTES
%   1. MATLAB's Rayleigh channel visualizes only up to (N/2+1)/fs in delay.
%   2. The actual "zero-time" in the CIR might be ambiguous.
%   3. Path attenuations are provided in dB but are plotted just as magnitude.
%
function [scenario, requiredBandwidth] = getMultipathScenario(scenarioType, scenarioIndex, fs)

  if nargin < 3 || isempty(fs)
    fs = 20e6; % Default sampling frequency is 20 MHz (802.11a/g)
  end

  % Test channel scenarios illustrating different multipath effects:
  % - Scenario 1: 1-tap channel. Flat magnitude and phase frequency response.
  % - Scenario 2: 2-tap channel. Adds fading and phase shift.
  % - Scenario 3: 2-tap channel. Second tap is one delay bin further compared to scenario 2.
  % - Scenario 4: 2-tap channel. Second tap is one delay bin further compared to scenario 3.
  % - Scenario 5: 4-tap channel.
  % - Scenario 6: 2-tap channel with aliased delay bin. between tap 2 and 3.
  % - Scenario 7: 3-tap channel.
  % - Scenario 8: 2-tap channel. Here, the second tap is at CP-1 length. The CFR will oscillate slowly.
  % - Scenario 9: 2-tap channel. Two closely spaced taps leading to merged taps.
  multipathScenariosTest = {
    struct('attenuations', 0, 'delays', 0/fs), ...
    struct('attenuations', [-3, -8], 'delays', [0, 1/fs]), ...
    struct('attenuations', [-3, -8], 'delays', [0, 2/fs]), ...
    struct('attenuations', [-3, -8], 'delays', [0, 3/fs]), ...
    struct('attenuations', [0, -3, -6, -9], 'delays', [0, 2/fs, 4/fs, 6/fs]), ...
    struct('attenuations', [-3, -8], 'delays', [0, 2.5/fs]), ...
    struct('attenuations', [-3, -8, -10], 'delays', [0, 0.25/fs, 1.5/fs]), ...
    struct('attenuations', [-3, -8], 'delays', [0, 15/fs]), ...
    struct('attenuations', [-3, -8], 'delays', [0, 1/50e6]), ...
  }; 

  % Realistic channel scenarios illustrating different static multipath effects for indoor Wi-Fi environments:
  % - Scenario 1: Simple direct path. Represents a line-of-sight (LoS) situation.
  % - Scenario 2: 2-tap channel. Mimics a direct path and a single reflection
  % - Scenario 3: 3-tap channel. Models a scenario with two significant reflections besides the LoS, such as from walls or furniture.
  % - Scenario 4: 4-tap channel. Simulates multiple reflections in a cluttered environment, like an office or living room.
  % - Scenario 5: 5-tap channel. Represents a complex indoor environment with multiple reflecting surfaces and objects.
  multipathScenariosReal = {
    struct('attenuations', [0], 'delays', [0]), ...
    struct('attenuations', [0, -6], 'delays', [0, 0.5/fs]), ...
    struct('attenuations', [0, -3, -9], 'delays', [0, 0.25/fs, 0.5/fs]), ...
    struct('attenuations', [0, -6, -12, -18], 'delays', [0, 0.15/fs, 0.3/fs, 0.45/fs]), ...
    struct('attenuations', [0, -3, -7, -10, -13, -15, -20, -25], 'delays', [0, 0.125/fs, 0.25/fs, 0.5/fs, 1/fs, 1.25/fs, 1.5/fs, 2/fs]), ...
  };


  if strcmp(scenarioType, 'test')
    multipathScenarios = multipathScenariosTest;
  else
    multipathScenarios = multipathScenariosReal;
  end

  scenario = multipathScenarios{scenarioIndex};
  fprintf('Selected scenario: %d\n', scenarioIndex);
  fprintf('\tAttenuations: %s\n', mat2str(scenario.attenuations));
  fprintf('\tDelays: %s\n', mat2str(scenario.delays));

  if any(scenario.delays > fs)
    error('Delays must be smaller than the sampling frequency!');
  end

  % Calculate bandwidth needed for resolution
  nonZeroDelays = scenario.delays(scenario.delays > 0); % Exclude zero delays as a 0 delay is always resolved
  if length(nonZeroDelays) > 0
    requiredBandwidth = computeRequiredBandwidth(scenario.delays, fs);
  else
    requiredBandwidth = fs;
  end

  fprintf('\tBandwidths needed for all paths: %.2f MHz\n', requiredBandwidth/1e6);
end


%computeRequiredBandwidth: Computes the combined bandwidth required for all delays in the given array.
%
% USAGE
%   requiredBandwidth = computeRequiredBandwidth(delays, fs);
%
% INPUT PARAMETERS
%   delays : Array of delays for which the bandwidth is to be computed.
%   fs     : Sampling frequency.
%
% OUTPUT PARAMETERS
%   requiredBandwidth : Bandwidth required to accommodate all paths.
%
% DETAILS
%   This function calculates the combined bandwidth required to accommodate all paths.
%
% REVISIT: Is this correct?
function requiredBandwidth = computeRequiredBandwidth(delays, fs)
  delaysInFs     = fs * delays;
  delaysInFsInt  = fix(delaysInFs);
  delaysInFsFrac = abs(delaysInFs - delaysInFsInt);
  nonZeroDelays  = delaysInFsFrac(delaysInFsFrac > 0); % These are the delays that can't be represented exactly with the given fs

  delayBandwidths = 1./(nonZeroDelays) * fs;

  if any(delaysInFsFrac == 0)
    delayBandwidths = [fs, delayBandwidths];
  end

  % Compute the least common multiple (LCM) of all bandwidths
  scSpacing    = 312.5e3; % We can only jump in multiples of 312.5 kHz
  lcmBandwidth = lcmm(floor(delayBandwidths * scSpacing));
  requiredBandwidth = ceil(lcmBandwidth / scSpacing);
end


function m = lcmm(args)
  % Compute LCM for a set of values (not limited to just two values)
  if isempty(args)
    m = 1;
    return
  end

  m = args(1);
  for i = 2:length(args)
    m = lcm(m, args(i));
  end
end
