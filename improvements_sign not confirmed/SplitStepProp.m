function BeamProp = SplitStepProp(Beam,KZ,distance)
%%  Propagate beam via FFT
% - "Beam" is the transverse structure
% - "KZ" is a grid of wave vector components in the propagation direction
%   (this is a k-space grid and is dependant on the wavelength)
%   (Note: This grid needs to be in the fftshifted form where the k-space 
%   origin is in the corners of the grid)
% - "distance" is the propagation distance we want to simulate the propagation over
% The output is the "Beam"'s transverse structure after a propagation "distance"



BeamFFT = fft2(Beam);
BeamK = BeamFFT.*exp(-1i*(KZ*distance));    
BeamProp = ifft2(BeamK);
