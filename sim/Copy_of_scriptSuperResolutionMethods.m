% scriptSuperResolutionMethods.m: This script is used to compare the MUSIC and MVDR super resolution methods.

clc;
clear;
close all;

%% Data loading
%------------------------------------------------------------------------------

dataFolder = '/Users/adamouazzani/Downloads/';
fileName  = 'label-path0_0cm_path1_690cm_2024-04-11_16-24-11_spliced.mat'; % 230 MHz

filePath = strcat(dataFolder, fileName);

load('test_splicing.mat');

whos;

%% Settings
%------------------------------------------------------------------------------

numPaths = 3;  % Number of paths to look for. We have 2 paths, but because the energy spreads itself a bit out, we may have to look further

mvdrEigEn     = true;  % If true, extract a subset of the data for MVDR  using singular value decomposition
mvdrEigCumSum = 0.997; % For simulated data

maxDistanceMeter = 15;      % How far to look for the signal
srDistResolution = 0.1;     % The desired distance resolution in meters
cableVelocity    = 0.7*3e8; % Speed of light in cable (approx)

%% Data Pre-processing
%------------------------------------------------------------------------------

csiArray = csiSpliced;
cirArray = cirSpliced;

activefreqList = infoSplice.ActiveFrequencies; % These are the frequencies of each sub-carrier/frequency in the channel frequency response, we need these later

nFrames = size(csiArray , 2);
n = size(csiArray, 1); % Number of frequency bins

% Reconstruct the CIR from the CSI using IFFT
reconstructedCIR = ifft(csiArray, [], 1);

% Plot the original and reconstructed CIR for the first frame (zoomed in)
figure;
numSamplesToShow = min(40, length(cirArray)); % Adjust this value as needed to zoom in on the peaks
stem(abs(cirArray(1:numSamplesToShow)), 'b', 'DisplayName', 'Original CIR', 'MarkerFaceColor', 'b');
hold on;
plot(abs(reconstructedCIR(1:numSamplesToShow)), 'r', 'DisplayName', 'Reconstructed CIR');
hold off;
title('Comparison of Original CIR and Reconstructed CIR from CSI');
xlabel('Sample Index');
ylabel('Magnitude');
legend;
grid on;

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
  csiSubset   = csiArray(:, i);
  csiForMusic = hankel(csiSubset);

  % 2) Covariance matrix estimation
  R = transpose(csiForMusic) * conj(csiForMusic) / size(csiForMusic , 1);
  Rinv = pinv(R);

  % 3) Perform eigendecomposition
  [eigVec, eigVal] = eig(R);
  eigVal = abs(diag(eigVal));

  % Ensure numPaths does not exceed the number of eigenvalues
  numPaths = min(numPaths, length(eigVal));

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
    mvdrRes = 1/(conj(transpose(steeringVec)) * RinvMvdr  * steeringVec);

    % MUSIC
    
    musicRes = 1 ./ sum(abs(transpose(conj(eigVecNoise)) * steeringVec) .^ 2);
    
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

% Combined across multiple measurements
cirArrayPwrAvg             = mean(cirArrayPwr, 1);
pseudoSpectrumMusicAvg     = mean(pseudoSpectrumMusic, 1);
pseudoSpectrumMvdrAvg      = mean(pseudoSpectrumMvdr, 1);
pseudoSpectrumMusicMvdrAvg = mean(pseudoSpectrumMusicMvdr,1);

% Ensure the averages have correct dimensions
cirArrayPwrAvg = cirArrayPwrAvg(:);  % Ensure it's a column vector
pseudoSpectrumMusicAvg = pseudoSpectrumMusicAvg(:);
pseudoSpectrumMvdrAvg = pseudoSpectrumMvdrAvg(:);
pseudoSpectrumMusicMvdrAvg = pseudoSpectrumMusicMvdrAvg(:);

cirArrayPwrAvg             = cirArrayPwrAvg ./ max(cirArrayPwrAvg, [], 'all');
pseudoSpectrumMusicAvg     = pseudoSpectrumMusicAvg ./ max(pseudoSpectrumMusicAvg, [], 'all');
pseudoSpectrumMvdrAvg      = pseudoSpectrumMvdrAvg ./ max(pseudoSpectrumMvdrAvg, [], 'all');
pseudoSpectrumMusicMvdrAvg = pseudoSpectrumMusicMvdrAvg ./ max(pseudoSpectrumMusicMvdrAvg, [], 'all');

%% Ensure cirArrayPwr has enough columns and at least 3 samples
numColsCirArrayPwr = size(cirArrayPwr, 2);
numColsToPlot = max(3, min(numColsCirArrayPwr, length(distanceAxis)));

% Ensure cirArrayPwrAvg has enough elements and at least 3 samples
numElemsCirArrayPwrAvg = length(cirArrayPwrAvg);
numElemsToPlot = max(3, min(numElemsCirArrayPwrAvg, length(distanceAxis)));

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

% Ensure cirArrayPwr has enough columns
numColsCirArrayPwr = size(cirArrayPwr, 2);
numColsToPlot = max(3, min(numColsCirArrayPwr, length(distanceAxis)));

% Plot the impulse response with distance
figure;
plot(distanceAxis(1:numColsToPlot), 10*log10(cirArrayPwr(1:numColsToPlot)));
xlabel('Distance (m)');
ylabel('Power (dB)');
title('Impulse response');
grid on;

% Plot the impulse response with distance for averages
figure;
hold on;
plot(distanceAxis(1:numElemsToPlot), 10*log10(cirArrayPwrAvg(1:numElemsToPlot)), 'LineWidth', 1.5, 'DisplayName', 'Impulse Response');
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
[cirPeaks, cirDistances] = findpeaks(10*log10(cirArrayPwrAvg(1:numElemsToPlot)), distanceAxis(1:numElemsToPlot), 'NPeaks', 2, 'SortStr', 'descend', 'MinPeakDistance', minPeakDistanceMeters);
if length(cirDistances) < 2
    warning('Not enough peaks found in the CIR for peak distance calculation.');
    cirPeakDistance = NaN;
else
    cirPeakDistance = cirDistances(2) - cirDistances(1);
end

% Find peaks and locations for pseudoSpectrumMusicAvg
[musicPeaks, musicDistances] = findpeaks(10*log10(pseudoSpectrumMusicAvg), srDistanceAxis, 'NPeaks', 2, 'SortStr', 'descend', 'MinPeakDistance', minPeakDistanceMeters);
if length(musicDistances) < 2
    warning('Not enough peaks found in the MUSIC pseudospectrum for peak distance calculation.');
    musicPeakDistance = NaN;
else
    musicPeakDistance = musicDistances(2) - musicDistances(1);
end

% Find peaks and locations for pseudoSpectrumMvdrAvg
[mvdrPeaks, mvdrDistances] = findpeaks(10*log10(pseudoSpectrumMvdrAvg), srDistanceAxis, 'NPeaks', 2, 'SortStr', 'descend', 'MinPeakDistance', minPeakDistanceMeters);
if length(mvdrDistances) < 2
    warning('Not enough peaks found in the MVDR pseudospectrum for peak distance calculation.');
    mvdrPeakDistance = NaN;
else
    mvdrPeakDistance = mvdrDistances(2) - mvdrDistances(1);
end

% Find peaks and locations for pseudoSpectrumMusicMvdrAvg
[musicMvdrPeaks, musicMvdrDistances] = findpeaks(10*log10(pseudoSpectrumMusicMvdrAvg), srDistanceAxis, 'NPeaks', 2, 'SortStr', 'descend', 'MinPeakDistance', minPeakDistanceMeters);
if length(musicMvdrDistances) < 2
    warning('Not enough peaks found in the MVDR + MUSIC pseudospectrum for peak distance calculation.');
    musicMvdrPeakDistance = NaN;
else
    musicMvdrPeakDistance = musicMvdrDistances(2) - musicMvdrDistances(1);
end

% Display the results
fprintf('Impulse Response:\n');
if ~isnan(cirPeakDistance)
    fprintf('Peak 1: %.2f dB at %.2f m\n', cirPeaks(1), cirDistances(1));
    fprintf('Peak 2: %.2f dB at %.2f m\n', cirPeaks(2), cirDistances(2));
    fprintf('Distance between peaks: %.2f m\n\n', cirPeakDistance);
else
    fprintf('Not enough peaks found for CIR.\n\n');
end

fprintf('MUSIC Pseudospectrum:\n');
if ~isnan(musicPeakDistance)
    fprintf('Peak 1: %.2f dB at %.2f m\n', musicPeaks(1), musicDistances(1));
    fprintf('Peak 2: %.2f dB at %.2f m\n', musicPeaks(2), musicDistances(2));
    fprintf('Distance between peaks: %.2f m\n\n', musicPeakDistance);
else
    fprintf('Not enough peaks found for MUSIC.\n\n');
end

fprintf('MVDR Pseudospectrum:\n');
if ~isnan(mvdrPeakDistance)
    fprintf('Peak 1: %.2f dB at %.2f m\n', mvdrPeaks(1), mvdrDistances(1));
    fprintf('Peak 2: %.2f dB at %.2f m\n', mvdrPeaks(2), mvdrDistances(2));
    fprintf('Distance between peaks: %.2f m\n\n', mvdrPeakDistance);
else
    fprintf('Not enough peaks found for MVDR.\n\n');
end

fprintf('MVDR + MUSIC Pseudospectrum:\n');
if ~isnan(musicMvdrPeakDistance)
    fprintf('Peak 1: %.2f dB at %.2f m\n', musicMvdrPeaks(1), musicMvdrDistances(1));
    fprintf('Peak 2: %.2f dB at %.2f m\n', musicMvdrPeaks(2), musicMvdrDistances(2));
    fprintf('Distance between peaks: %.2f m\n', musicMvdrPeakDistance);
else
    fprintf('Not enough peaks found for MVDR + MUSIC.\n');
end

%% Peak extraction and delay estimation variance
%------------------------------------------------------------------------------

numMeasurements = length(cirArray); % cirArray is a vector
delayEstimates = zeros(numMeasurements, 3); % [CIR, MUSIC, MVDR]

minPeakDistanceMeters = 2; % Minimum peak distance in meters

for i = 1:numMeasurements
    % CIR Peaks
    [cirPeaks, cirDistances] = findpeaks(10*log10(abs(cirArray(1:length(distanceAxis)))), distanceAxis, ...
                                         'NPeaks', numPaths, 'SortStr', 'descend', 'MinPeakDistance', minPeakDistanceMeters);
    if length(cirDistances) == 2
        delayEstimates(i, 1) = cirDistances(2) - cirDistances(1);
    end

    % MUSIC Peaks
    if size(pseudoSpectrumMusic, 1) >= i
        [musicPeaks, musicDistances] = findpeaks(10*log10(pseudoSpectrumMusic(i, :)), srDistanceAxis, ...
                                                 'NPeaks', numPaths, 'SortStr', 'descend', 'MinPeakDistance', minPeakDistanceMeters);
        if length(musicDistances) == 2
            delayEstimates(i, 2) = musicDistances(2) - musicDistances(1);
        end
    end

    % MVDR Peaks
    if size(pseudoSpectrumMvdr, 1) >= i
        [mvdrPeaks, mvdrDistances] = findpeaks(10*log10(pseudoSpectrumMvdr(i, :)), srDistanceAxis, ...
                                               'NPeaks', numPaths, 'SortStr', 'descend', 'MinPeakDistance', minPeakDistanceMeters);
        if length(mvdrDistances) == 2
            delayEstimates(i, 3) = mvdrDistances(2) - mvdrDistances(1);
        end
    end
end

% Remove NaN values from delay estimates
delayEstimates = delayEstimates(~any(isnan(delayEstimates), 2), :);

% Boxplot of delay estimation variance
figure;
boxplot(delayEstimates, 'Labels', {'CIR', 'MUSIC', 'MVDR'});
ylabel('Delay Estimate (m)');
title('Variance of Delay Estimation between Two Paths');
grid on;

% Estimate true delay from original CIR (using the first frame as an example)
[cirPeaks, cirLocations] = findpeaks(abs(cirArray), 'NPeaks', 2, 'SortStr', 'descend', 'MinPeakDistance', 2);
trueDelay = abs(cirLocations(2) - cirLocations(1));

fprintf('Estimated True Delay: %.2f samples\n', trueDelay);

% Number of runs for variance and MSE calculation
nRuns = 1;

% Initialize arrays to store delay estimates
cirPeakDistances = zeros(nRuns, 1);
musicPeakDistances = zeros(nRuns, 1);
mvdrPeakDistances = zeros(nRuns, 1);

% Initialize arrays to store MSE values
mseCIR = zeros(nRuns, 1);
mseMusic = zeros(nRuns, 1);
mseMvdr = zeros(nRuns, 1);

% Loop over runs
for i = 1:nRuns
    
    noisyCSI = csiArray + 0.00 * randn(size(csiArray));

    % Reconstruct CIR from noisy CSI
    reconstructedCIR = ifft(noisyCSI, [], 1);

    % Find peaks for CIR
    [cirPeaks, cirLocations] = findpeaks(abs(reconstructedCIR(:, 1)), 'NPeaks', 2, 'SortStr', 'descend', 'MinPeakDistance', 2);
    cirPeakDistances(i) = abs(cirLocations(2) - cirLocations(1));

    % Calculate MUSIC pseudospectrum for the current noisy CSI
    pseudoSpectrumMusicCurrent = zeros(1, numel(srTimeDelayAxis));
    pseudoSpectrumMvdrCurrent = zeros(1, numel(srTimeDelayAxis));

    % Form the Hankel matrix
    csiSubset = noisyCSI(:, 1);
    csiForMusic = hankel(csiSubset);

    % Covariance matrix estimation
    R = transpose(csiForMusic) * conj(csiForMusic) / size(csiForMusic, 1);
    Rinv = pinv(R);

    % Perform eigendecomposition
    [eigVec, eigVal] = eig(R);
    eigVal = abs(diag(eigVal));

    % Get the signal and noise subspaces
    [eigValSort, idx] = sort(eigVal, 'descend');
    eigVecNoise = eigVec(:, idx(numPaths + 1:end));

    if mvdrEigEn
        % Find the indices of the eigenvalues that are below the threshold
        eigSelectIndices = idx(cumsum(eigValSort / sum(eigValSort)) < mvdrEigCumSum);

        % Transform the CSI data to the signal subspace
        signalSubspace = eigVec(:, eigSelectIndices);
        transformedCSI = csiForMusic * signalSubspace;
        transformedCSI = transpose(transformedCSI);

        % MVDR
        Rmvdr = transpose(transformedCSI) * conj(transformedCSI) / size(transformedCSI, 1);
        RinvMvdr = pinv(Rmvdr);
    else
        RinvMvdr = Rinv;
    end

    % Compute the pseudospectrum
    for j = 1:numel(srTimeDelayAxis)
        steeringVec = exp(-1i * 2 * pi * srTimeDelayAxis(j) * activefreqList);

        % MVDR from (spatial filtering)
        mvdrRes = 1 / (conj(transpose(steeringVec)) * RinvMvdr * steeringVec);

        % MUSIC
        musicRes = 1 / sum(abs(transpose(conj(eigVecNoise)) * steeringVec).^2);

        pseudoSpectrumMvdrCurrent(j) = abs(mvdrRes)^2;
        pseudoSpectrumMusicCurrent(j) = musicRes;
    end

    % Find peaks for MUSIC
    [musicPeaks, musicLocations] = findpeaks(10*log10(pseudoSpectrumMusicCurrent), srDistanceAxis, 'NPeaks', 2, 'SortStr', 'descend', 'MinPeakDistance', minPeakDistanceMeters);
    musicPeakDistances(i) = abs(musicLocations(2) - musicLocations(1));

    % Find peaks for MVDR
    [mvdrPeaks, mvdrLocations] = findpeaks(10*log10(pseudoSpectrumMvdrCurrent), srDistanceAxis, 'NPeaks', 2, 'SortStr', 'descend', 'MinPeakDistance', minPeakDistanceMeters);
    mvdrPeakDistances(i) = abs(mvdrLocations(2) - mvdrLocations(1));

    % Calculate MSE
    mseCIR(i) = (cirPeakDistances(i) - trueDelay)^2;
    mseMusic(i) = (musicPeakDistances(i) - trueDelay)^2;
    mseMvdr(i) = (mvdrPeakDistances(i) - trueDelay)^2;
end

% Plot the MSE results
figure;
bar([mean(mseCIR), mean(mseMusic), mean(mseMvdr)]);
xticklabels({'CIR', 'MUSIC', 'MVDR'});
xlabel('Method');
ylabel('Mean Square Error (m^2)');
title('Mean Square Error of Delay Estimation');
grid on;
