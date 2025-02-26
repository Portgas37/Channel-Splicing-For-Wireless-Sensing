%getPrunedDFTMatrix: Generates a pruned DFT matrix where only the active FFT indices are retained.
% This allows you to work with systems that do not use all subcarriers and estimate the CIR or CFR
%
% USAGE
%   prunedDFTMat = getPrunedDFTMatrix(nfft, activeFFTIndices, fftDcIdx)
%
% INPUT PARAMETERS
%   nfft             : Length of the FFT
%   activeFFTIndices : Indices of active FFT subcarriers to be pruned. This assumes that the DC is in the center
%   fftDcIdx         : Index for the carrier in the FFT, typically used to ensure correct DC positioning.
%
% OUTPUT PARAMETERS
%   prunedDFTMat: The pruned DFT matrix based on the active FFT indices. For multiplication, it assumes the frequency domain
%                 signal has the same length as the active FFT indices and that it matches these indices.
%                 The output product will be the impulse response starting at time 0.
%
% DETAILS
%   This function produces a pruned DFT matrix based on specified active FFT indices.
%   The main intention behind pruning the DFT matrix is to work with systems that do not use all subcarriers.
%
%   Typically, the FFT size or corresponds to factors like channel bandwidth and subcarrier spacing.
%   For a 20 MHz channel bandwidth with a 312.5 KHz subcarrier spacing, the FFT size can be 64. But, if
%   there is an additional subcarrier, the size might increase to 65 or you can keep it to size 64
%   and not use the last sub-carrier, you just have to keep the sub-carrier spacing fixed.
%
function [prunedDFTMat] = getPrunedDFTMatrix(nfft, activeFFTIndices, fftDcIdx)

  if nargin < 3
    fftDcIdx = floor(nfft / 2) + 1;
  end

  % Shift the DFT matrix rows to center the DC component along the rows (in columns, DC is still at the first column)
  dftMatCenter = fftshift(dftmtx(nfft), 1);
  prunedDFTMat = dftMatCenter(activeFFTIndices, :);

  % Ensure the centered DFT matrix's carrier index is all ones for DC integrity
  if ~all(dftMatCenter(fftDcIdx, :) == 1)
    error('The FFT carrier index row is expected to contain all ones for correct DC positioning.');
  end
end
