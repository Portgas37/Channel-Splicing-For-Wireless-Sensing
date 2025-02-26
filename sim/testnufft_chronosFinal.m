%testnufft.m
%
% This script tests use of the non-uniform FFT functions
%
% References:
% - https://ch.mathworks.com/matlabcentral/answers/2032659-adjoint-inverse-nufft

clear; clc;
close all;
rng("default");
3
load('test_splicing.mat'); % Load the data from a previous run of scriptCsiMultiBand.m
                           % When you run scriptCsiMultiBand.m, a file called test_splicing.mat is created

                       

%% Test the data
%==============================================================================
fignum = 22; % Some large figure number to not overwrite previous figures if you have any other plots open

% For the frequency we use the baseband frequencies, for the time, we use a subset of the CIR
freqAxis = infoSplice.ActiveFrequencies;
timeAxis = (0:(infoSplice.FFTLength)-1) / (infoSplice.SampleRate);
timeAxis = timeAxis(timeAxis < 8/20e6); % Only keep up to 8/20e6=0.4us in time, so 8 main taps when in 20 MHz
% Print freqAxis

% TODO: Add gaps here

% We can get the CIR from the CSI data by using nufft with the frequency and time axis (negative)
% When going from frequency to time with nufft, you have to divide by the length of the input to nufft to get the correct scaling
% This is not the case when going from time to frequency, here, there is no scaling needed
cirRef = nufft(csiSpliced, freqAxis, -timeAxis) / length(csiSpliced);

figure(fignum); 
stem(timeAxis, abs(cirRef));
fignum = fignum + 1;
title('Channel impulse response reconstructed using NUFFT')

%% Non-uniform DFT matrix
%==============================================================================

% Construct the non-uniform DFT matrix F
% Note that in the "Decimeter-Level Localization with a Single WiFi Access Point" paper
% they need ||F||^2, we might not have to calculate this. Since each element in the Fourier matrix
% is a complex exponential, it will have magnitude 1, so ||F||_2 is a scaling
% to correct for the size of the changing when doing the transformation
F = exp(-1j*2*pi*reshape(timeAxis, [], 1)*reshape(freqAxis, 1, []));
Finv = pinv(F);
norm(F,2) % ||F||_2, this is just 16

% Verify that you get a smoother version but approximately correct version of the CSI
% Note that when we created the CIR, we trimmed it to 0.2us, so we only have 4 main taps
% Because these are the main taps, and most later taps are noise, using DFT here also filters the nosie
csiFiltered = transpose(F) * cirRef;
figure(fignum);
plot(abs(csiFiltered))
hold on;
plot(abs(csiSpliced))
hold off;
fignum = fignum + 1;
title('Comparison of original CSI with CSI from the trimmed CIR')

cirFiltered = transpose(Finv) * csiFiltered;
figure(fignum);
stem(abs(cirFiltered))
hold on;
plot(abs(cirRef))
hold off;
fignum = fignum + 1;
title('Check that we can use the inverted non-uniform DFT matrix to get the CIR from the csi')


alpha = 5; % Adjust alpha value
epsilon = 0.001; % Adjust epsilon value

% Performing the channel cut
%N = length(csiSpliced);
%idxCut = fix(N/2);
%idxCut = 1;
%csiSplicedCut = csiSpliced;
%csiSplicedCut(idxCut:idxCut+810) = [];
%freqAxisCut = freqAxis;
%freqAxisCut(idxCut:idxCut+810) = [];
%Fcut = exp(-1j*2*pi*reshape(timeAxis, [], 1)*reshape(freqAxisCut, 1, []));

% Initialize variables
N = length(csiSpliced);
csiSplicedCut = csiSpliced;
freqAxisCut = freqAxis;

% Remove the first 600 and the last 600 elements
keepIndices = true(N, 1);
keepIndices(1:600) = false;
keepIndices(end-599:end) = false;

% Define the number of elements to remove and keep for the remaining part
removeCount = 100;
keepCount = 200;

% Adjust the starting index after the initial removal
startIndex = 600;
endIndex = N - 600;

% Loop through the data to set indices for uniform removal
for i = startIndex:(removeCount + keepCount):endIndex
    % Calculate the end index for removal within the allowed range
    uniformEndIndex = min(i + removeCount - 1, endIndex);
    
    % Set the indices to false for the elements to remove
    keepIndices(i:uniformEndIndex) = false;
end

% Apply the keepIndices to both csiSpliced and freqAxis
csiSplicedCut = csiSpliced(keepIndices);
freqAxisCut = freqAxis(keepIndices);

% Construct the Fcut matrix with the updated freqAxisCut
Fcut = exp(-1j*2*pi*reshape(timeAxis, [], 1)*reshape(freqAxisCut, 1, []));



% Choosing the number of runs
N_runs = 30;

% Initialize storage for magnitudes and phases of reconstructed p
reconstructed_magnitudes = zeros(N_runs, length(cirRef)); 
reconstructed_phases = zeros(N_runs, length(cirRef));

% Loop for computation and storing results
for runIndex = 1:N_runs
    rng(runIndex); 

    p = computeInverseNDFT(csiSplicedCut, Fcut, alpha, epsilon);

    % Store magnitude and phase of reconstructed p
    reconstructed_magnitudes(runIndex, :) = abs(p);
    reconstructed_phases(runIndex, :) = angle(p);
end

% Average reconstructed p over all runs
average_reconstructed_p = mean(reconstructed_magnitudes .* exp(1i * reconstructed_phases), 1);

% Plot the average reconstructed p
figure;
stem(timeAxis, abs(average_reconstructed_p), 'b', 'DisplayName', 'Reconstructed CIR from p');
hold on;
plot(timeAxis, abs(cirRef), 'r', 'DisplayName', 'Original CIR');
title('Comparison of Original CIR and Reconstructed CIR from p');
xlabel('Time (s)');
ylabel('Magnitude');
legend;
grid on;

% Define the speed of light in meters per second
c = 3e8; % speed of light in m/s

% Convert time axis to distance axis
distanceAxis = timeAxis * c / 2; % divide by 2 for round-trip distance

% Plot the CIR with distance axis
figure;
stem(distanceAxis, abs(average_reconstructed_p), 'b', 'DisplayName', 'Reconstructed CIR from p');
hold on;
plot(distanceAxis, abs(cirRef), 'r', 'DisplayName', 'Original CIR');
title('Comparison of Original CIR and Reconstructed CIR from p');
xlabel('Distance (m)');
ylabel('Magnitude');
legend;
grid on;

% Plotting all magnitudes and phases as points with distance axis
figure; 
subplot(2,1,1); plot(distanceAxis, reconstructed_magnitudes', '.', 'MarkerSize', 12); title('Reconstructed Magnitudes for All Runs');
subplot(2,1,2); plot(distanceAxis, reconstructed_phases', '.', 'MarkerSize', 12); title('Reconstructed Phases for All Runs');


% Plotting all magnitudes and phases as points
figure; 
subplot(2,1,1); plot(reconstructed_magnitudes', '.', 'MarkerSize', 12); title('Reconstructed Magnitudes for All Runs');
subplot(2,1,2); plot(reconstructed_phases', '.', 'MarkerSize', 12); title('Reconstructed Phases for All Runs');

% True magnitude and phase for comparison
true_magnitude = abs(cirRef);
true_phase = angle(cirRef);

% Initialize loss arrays
magnitude_losses = zeros(N_runs, 1);
phase_losses = zeros(N_runs, 1);

% Compute squared relative losses for magnitude and phase
for i = 1:N_runs
    % Magnitude Loss Calculation
    squared_diffs_mag = (reconstructed_magnitudes(i, :) - true_magnitude').^2; % Ensure dimension match
    squared_relative_diffs_mag = squared_diffs_mag ./ (true_magnitude').^2;
    magnitude_losses(i) = sum(squared_relative_diffs_mag) / numel(true_magnitude);
    
    % Phase Loss Calculation
    % Normalize phase differences within [-pi, pi] range
    phase_diff = wrapToPi(reconstructed_phases(i, :) - true_phase');
    squared_diffs_phase = phase_diff.^2; % Squared differences
    squared_relative_diffs_phase = squared_diffs_phase ./ (true_phase'.^2);
    phase_losses(i) = sum(squared_relative_diffs_phase) / length(true_phase);
end

% Compute mean loss for all runs
mean_magnitude_loss = mean(magnitude_losses);
mean_phase_loss = mean(phase_losses);
fprintf('Mean Magnitude Loss: %f\n', mean_magnitude_loss);
fprintf('Mean Phase Loss: %f\n', mean_phase_loss);

% Example: Finding the frequency corresponding to the 600th index
index = 600;
if index <= length(freqAxis)
    frequencyValue = freqAxis(index);
    fprintf('The frequency corresponding to index %d is %.2f MHz\n', index, frequencyValue / 1e6);
else
    fprintf('Index %d is out of bounds for the frequency axis.\n', index);
end



function p = computeInverseNDFT(h, F, alpha, epsilon)
    p = randn(size(F, 1), 1);
    t = 0;

    gamma = 0.0001;
    converged = false;

    while ~converged
        p_next = p - gamma * (F').' * ((F.' * p) - h);
        p_next = sparsify(p_next, gamma * alpha);

        fprintf('Iteration: %d, Norm difference: %f\n', t, norm(p_next - p));

        if norm(p_next - p) < epsilon
            converged = true;
        end

        p = p_next;
        t = t + 1;

        % Add a check for convergence to avoid infinite loop
        if t > 400 % Maximum number of iterations
            fprintf('Maximum number of iterations reached.\n');
            break;
        end
        
    end
end

function p = sparsify(p, threshold)
    for i = 1:length(p)
        if abs(p(i)) < threshold
            p(i) = 0;
        else
            p(i) = p(i) * (abs(p(i)) - threshold) / abs(p(i));
        end
    end
end