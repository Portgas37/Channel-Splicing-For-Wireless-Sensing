%generateFreqAxis: Generates a frequency axis for given sampling frequency, FFT size, and DC index.
% The frequency axis starts from the most negative frequency and goes to most positive with the
% DC component at the center.
%
% INPUTS:
%   fs   : Sampling frequency (Hz)
%   nfft : FFT size (number of points)
%
% OUTPUT:
%   freqVec : Frequency vector in MHz
%
function freqVec = generateFreqAxis(fs, nfft)
  fftDcIdx = getFftDcIdx(nfft);
  freqVec  = fs / nfft * (-(fftDcIdx-1) : nfft - fftDcIdx);
end

function fftDcIdx = getFftDcIdx(nfft)
  fftDcIdx = floor(nfft / 2) + 1;
end
