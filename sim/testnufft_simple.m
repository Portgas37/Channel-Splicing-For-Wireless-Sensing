%testnufft.m
%
% This script tests use of the non-uniform FFT functions
%
% References:
% - https://ch.mathworks.com/matlabcentral/answers/2032659-adjoint-inverse-nufft

clear; clc;
close all;
rng("default");

load('test_splicing.mat'); % Load the data from a previous run of scriptCsiMultiBand.m
                           % When you run scriptCsiMultiBand.m, a file called test_splicing.mat is created

%% Test the data
%==============================================================================
fignum = 22; % Some large figure number to not overwrite previous figures if you have any other plots open

% For the frequency we use the baseband frequencies, for the time, we use a subset of the CIR
freqAxis = infoSplice.ActiveFrequencies;
timeAxis = (0:(infoSplice.FFTLength)-1) / (infoSplice.SampleRate);
timeAxis = timeAxis(timeAxis < 8/20e6); % Only keep up to 8/20e6=0.4us in time, so 8 main taps when in 20 MHz

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
