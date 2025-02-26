%spliceSettingsCommon.m: Common settings for the Wi-Fi transmitter and receiver simulation setup

% WLAN configuration (we use 802.11n, which is what the HT indicates)
cfgHT = wlanHTConfig( ...
  'ChannelBandwidth', 'CBW20', ...
  'NumTransmitAntennas', 1, ...
  'RecommendSmoothing', false, ...
  'PSDULength', 0);

ind = wlanFieldIndices(cfgHT); % Wi-Fi frame field indices for extracting L-LTF and HT-LTF fields
fs  = wlanSampleRate(cfgHT);

infoLLTF  = wlanHTOFDMInfo('L-LTF', cfgHT); % L-LTF information. infoLLTF describes which sub-carriers are active etc.
scSpacing = fs / infoLLTF.FFTLength;
nfft      = infoLLTF.FFTLength;

infoHTLTF = wlanHTOFDMInfo('HT-LTF', cfgHT); % HT-LTF information
