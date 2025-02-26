% scriptSuperResolutionMethods.m: This script is used to compare the MUSIC and MVDR super resolution methods.

clc;
clear;
close all;

%% Data loading
%------------------------------------------------------------------------------

% TODO: Pick your dataset and make the right path on your computer
dataFolder = '/Users/adamouazzani/Downloads/';
%fileName  = 'label-path0_0cm_path1_690cm_2024-04-11_16-24-11_spliced.mat'; % 230 MHz
fileName  = 'label-path0_0cm_path1_690cm_2024-04-11_16-51-20_spliced.mat';  % 80 MHz


filePath = strcat(dataFolder, fileName);

load(filePath);

%% Settings
%------------------------------------------------------------------------------

% TODO: Try to understand these parameters
numPaths = 3;  % Number of paths to look for. We have 2 paths, but because the energy spreads itself a bit out, we may have to look further

mvdrEigEn     = true;  % If true, extract a subset of the data for MVDR  using singular value decomposition
%mvdrEigCumSum = 0.998; % Around 0.998 seems stable for large 230 MHz BW
mvdrEigCumSum = 0.995; % Around 0.995 seems stable for smaller 80 MHz BW

maxDistanceMeter = 15;      % How far to look for the signal
srDistResolution = 0.1;     % The desired distance resolution in meters
cableVelocity    = 0.7*3e8; % Speed of light in cable (approx)

%% Data Pre-processing
%------------------------------------------------------------------------------

% TODO: Here, you can see how I get the data out
csiArray   = dataStruct.csiSplicedArray;
cirArray   = dataStruct.cirSplicedArray;
infoSplice = dataStruct.infoSplice;

activefreqList = infoSplice.ActiveFrequencies; % These are the frequencies of each sub-carrier/frequency in the channel frequency response, we need these later

nFrames = 1; % TODO: Number of frames to process, you can put 1 or the number of frames you have
%nFrames = size(csiArray , 1);


%% Super resolution MUSIC using multiple snapshot
%------------------------------------------------------------------------------


distResolution       = cableVelocity/infoSplice.SampleRate;
srTimeResolutionGain = (distResolution/srDistResolution); % How much to increase the time resolution from just the sample rate for the super-resolution methods
srTimeResolution     = 1/(srTimeResolutionGain*infoSplice.SampleRate); % Time resolution for super-resolution methods
srDistResolution     = srTimeResolution * cableVelocity;               % Distance resolution for super-resolution methods
maxTimeDelay         = maxDistanceMeter/cableVelocity;

distanceAxis    = 0:cableVelocity/infoSplice.SampleRate:maxDistanceMeter; % Distance axis for just the impulse response
srDistanceAxis  = 0:srDistResolution:(maxDistanceMeter+srDistResolution); % Super-resolution distance axis
srTimeDelayAxis = srDistanceAxis / cableVelocity;

pseudoSpectrumMusic = zeros(nFrames, numel(srTimeDelayAxis), 1); % Pseudospectrum for MUSIC
pseudoSpectrumMvdr  = zeros(nFrames, numel(srTimeDelayAxis), 1); % Pseudospectrum for MVDR

for i = 1:nFrames

  % 1) Form the Hankel matrix
  csiSubset   = csiArray(i, :);
  csiForMusic = hankel(csiSubset);

  % 2) Covariance matrix estimation
  R = transpose(csiForMusic) * conj(csiForMusic) / size(csiForMusic , 1);
  Rinv = pinv(R);

  % 3) Perform eigendecomposition
  [eigVec, eigVal] = eig(R);
  eigVal = abs(diag(eigVal));

  % 4) Get the signal and noise subspaces
  [eigValSort, idx] = sort(eigVal, 'descend');
  eigVecSig   = eigVec(:, idx(1:numPaths));
  eigVecNoise = eigVec(:, idx(numPaths+1:end));

  if mvdrEigEn

    % Find the indices of the eigenvalues that are below the threshold, this allows us to select a subset of the data that has the most energy
    % while trying to get rid of noise
    eigSelectIndices = idx(cumsum(eigValSort/sum(eigValSort)) < mvdrEigCumSum);

    % Transform the CSI data to the signal subspace
    signalSubspace = eigVec(:, eigSelectIndices);
    transformedCSI = csiForMusic * signalSubspace;
    transformedCSI = transpose(transformedCSI);

    % MVDR in general has the peak more in the correct place, but is less steep
    Rmvdr    = transpose(transformedCSI) * conj(transformedCSI) / size(transformedCSI , 1);
    RinvMvdr = pinv(Rmvdr);

  else
    RinvMvdr = Rinv;
  end

  % 5) Compute the pseudospectrum
  for j = 1:numel(srTimeDelayAxis)
    steeringVec = exp(-1i * 2 * pi * srTimeDelayAxis(j) * activefreqList);

    % MVDR from (spatial filtering)
    % See the paper Multi-Band Superresolution Multipath Channel Path Delay Estimation for CIR-Based Localization
    mvdrRes = 1/(conj(transpose(steeringVec)) * RinvMvdr  * steeringVec);

    % MUSIC
    musicRes = 1./sum(abs(transpose(conj(eigVecNoise)) * steeringVec).^2);

    pseudoSpectrumMvdr(i, j)  = abs(mvdrRes)^2;
    pseudoSpectrumMusic(i, j) = musicRes;
  end
end

cirArrayPwr             = abs(cirArray).^2;
pseudoSpectrumMusicMvdr = pseudoSpectrumMvdr .* pseudoSpectrumMusic;

% We divide by the max to always have the strongest component at 0 dB
cirArrayPwr             = cirArrayPwr ./ max(cirArrayPwr, [], 'all');
pseudoSpectrumMvdr      = pseudoSpectrumMvdr ./ max(pseudoSpectrumMvdr, [], 'all');
pseudoSpectrumMusic     = pseudoSpectrumMusic ./ max(pseudoSpectrumMusic, [], 'all');
pseudoSpectrumMusicMvdr = pseudoSpectrumMusicMvdr ./ max(pseudoSpectrumMusicMvdr, [], 'all');

%% Plots

% Plot the pseudospectrum with distance
figure;
plot(srDistanceAxis, 10*log10(pseudoSpectrumMusic));
xlabel('Distance (m)');
ylabel('Pseudospectrum (dB)');
title('MUSIC Pseudospectrum');
grid on

% Plot the pseudospectrum with distance
figure;
plot(srDistanceAxis, 10*log10(pseudoSpectrumMvdr));
xlabel('Distance (m)');
ylabel('Pseudospectrum (dB)');
title('MVDR Pseudospectrum');
grid on

% Plot the pseudospectrum with distance
figure;
plot(srDistanceAxis, 10*log10(pseudoSpectrumMusicMvdr));
xlabel('Distance (m)');
ylabel('Pseudospectrum (dB)');
title('MVDR + MUSIC Pseudospectrum');
grid on

% Plot the impulse response with distance
figure;
plot(distanceAxis, 10*log10(cirArrayPwr(:,1:length(distanceAxis))));
xlabel('Distance (m)');
ylabel('Power (dB)');
title('Impulse response');
grid on

% Combined across multiple measurements
cirArrayPwrAvg             = mean(cirArrayPwr, 1);
pseudoSpectrumMusicAvg     = mean(pseudoSpectrumMusic, 1);
pseudoSpectrumMvdrAvg      = mean(pseudoSpectrumMvdr, 1);
pseudoSpectrumMusicMvdrAvg = mean(pseudoSpectrumMusicMvdr,1);

cirArrayPwrAvg             = cirArrayPwrAvg ./ max(cirArrayPwrAvg, [], 'all');
pseudoSpectrumMusicAvg     = pseudoSpectrumMusicAvg ./ max(pseudoSpectrumMusicAvg, [], 'all');
pseudoSpectrumMvdrAvg      = pseudoSpectrumMvdrAvg ./ max(pseudoSpectrumMvdrAvg, [], 'all');
pseudoSpectrumMusicMvdrAvg = pseudoSpectrumMusicMvdrAvg ./ max(pseudoSpectrumMusicMvdrAvg, [], 'all');

figure;
hold on;
plot(distanceAxis, 10*log10(cirArrayPwrAvg(1:length(distanceAxis))), 'LineWidth', 1.5, 'DisplayName', 'Impulse Response');
plot(srDistanceAxis, 10*log10(pseudoSpectrumMusicAvg), 'LineWidth', 1.5, 'DisplayName', 'MUSIC Pseudospectrum');
plot(srDistanceAxis, 10*log10(pseudoSpectrumMvdrAvg), 'LineWidth', 1.5, 'DisplayName', 'MVDR Pseudospectrum');
plot(srDistanceAxis, 10*log10(pseudoSpectrumMusicMvdrAvg), 'LineWidth', 1.5, 'DisplayName', 'MVDR + MUSIC Pseudospectrum');
xlabel('Distance (m)');
ylabel('Power (dB)');
title({'Comparison of Impulse Response and Pseudospectra'; 'for Super-Resolution Delay Estimation'});
legend('Location', 'best');
grid on;

%% Peak extraction

% Settings
%------------------------------------------------------------------------------
minPeakDistanceMeters = 2; % Minimum peak distance in meters

% Find peaks and locations for cirArrayPwrAvg
[cirPeaks, cirDistances] = findpeaks(10*log10(cirArrayPwrAvg(1:length(distanceAxis))), distanceAxis, 'NPeaks', 2, 'SortStr', 'descend', 'MinPeakDistance', minPeakDistanceMeters);
cirPeakDistance = cirDistances(2) - cirDistances(1);

% Find peaks and locations for pseudoSpectrumMusicAvg
[musicPeaks, musicDistances] = findpeaks(10*log10(pseudoSpectrumMusicAvg), srDistanceAxis, 'NPeaks', 2, 'SortStr', 'descend', 'MinPeakDistance', minPeakDistanceMeters);
musicPeakDistance = musicDistances(2) - musicDistances(1);

% Find peaks and locations for pseudoSpectrumMvdrAvg
[mvdrPeaks, mvdrDistances] = findpeaks(10*log10(pseudoSpectrumMvdrAvg), srDistanceAxis, 'NPeaks', 2, 'SortStr', 'descend', 'MinPeakDistance', minPeakDistanceMeters);
mvdrPeakDistance = mvdrDistances(2) - mvdrDistances(1);

% Find peaks and locations for pseudoSpectrumMusicMvdrAvg
[musicMvdrPeaks, musicMvdrDistances] = findpeaks(10*log10(pseudoSpectrumMusicMvdrAvg), srDistanceAxis, 'NPeaks', 2, 'SortStr', 'descend', 'MinPeakDistance', minPeakDistanceMeters);
musicMvdrPeakDistance = musicMvdrDistances(2) - musicMvdrDistances(1);

% Display the results
fprintf('Impulse Response:\n');
fprintf('Peak 1: %.2f dB at %.2f m\n', cirPeaks(1), cirDistances(1));
fprintf('Peak 2: %.2f dB at %.2f m\n', cirPeaks(2), cirDistances(2));
fprintf('Distance between peaks: %.2f m\n\n', cirPeakDistance);

fprintf('MUSIC Pseudospectrum:\n');
fprintf('Peak 1: %.2f dB at %.2f m\n', musicPeaks(1), musicDistances(1));
fprintf('Peak 2: %.2f dB at %.2f m\n', musicPeaks(2), musicDistances(2));
fprintf('Distance between peaks: %.2f m\n\n', musicPeakDistance);

fprintf('MVDR Pseudospectrum:\n');
fprintf('Peak 1: %.2f dB at %.2f m\n', mvdrPeaks(1), mvdrDistances(1));
fprintf('Peak 2: %.2f dB at %.2f m\n', mvdrPeaks(2), mvdrDistances(2));
fprintf('Distance between peaks: %.2f m\n\n', mvdrPeakDistance);

fprintf('MVDR + MUSIC Pseudospectrum:\n');
fprintf('Peak 1: %.2f dB at %.2f m\n', musicMvdrPeaks(1), musicMvdrDistances(1));
fprintf('Peak 2: %.2f dB at %.2f m\n', musicMvdrPeaks(2), musicMvdrDistances(2));
fprintf('Distance between peaks: %.2f m\n', musicMvdrPeakDistance);
