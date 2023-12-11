%function BeamProp = SplitStepProp(Beam,KZ,distance)
%%  Propagate beam via FFT
% - "Beam" is the transverse structure
% - "KZ" is a grid of wave vector components in the propagation direction
%   (this is a k-space grid and is dependant on the wavelength)
%   (Note: This grid needs to be in the fftshifted form where the k-space 
%   origin is in the corners of the grid)
% - "distance" is the propagation distance we want to simulate the propagation over
% The output is the "Beam"'s transverse structure after a propagation "distance"
%{
L1= 5e-3 ;
xinstart=-L1./2;
xinend=L1./2;
yinstart=-L1./2;
yinend=L1./2; 
mxin=5e-3./0.008e-3;
myin=mxin;
xstart=xinstart;
xend=xinend;
ystart=yinstart;
yend=yinstart;
mxout=mxin;
myout=myin;

[BeamProp,pixelout] = Scalar_Bluestein(Beam,mxin,myin,0.008e-3,distance,xstart,xend,ystart,yend,mxout,myout);
BeamProp=BeamProp;
%}
