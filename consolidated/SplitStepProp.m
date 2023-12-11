function BeamProp = SplitStepProp(Beam,distance)
%%  Propagate beam via FFT
% - "Beam" is the transverse structure
% - "KZ" is a grid of wave vector components in the propagation direction
%   (this is a k-space grid and is dependant on the wavelength)
%   (Note: This grid needs to be in the fftshifted form where the k-space 
%   origin is in the corners of the grid)
% - "distance" is the propagation distance we want to simulate the propagation over
% The output is the "Beam"'s transverse structure after a propagation "distance"
 
global WidthX WidthY

xinstart=-WidthX./2;
xinend=WidthX./2;
yinstart=-WidthY./2;
yinend=WidthY./2; 

global PixelSize

mxin=(WidthX)./(PixelSize);
myin=mxin;
xstart=xinstart;
xend=xinend;
ystart=yinstart;
yend=yinstart;
mxout=mxin;
myout=myin;

[BeamProp,pixelout] = Scalar_Bluestein(Beam,mxin,myin,PixelSize,distance,xstart,xend,ystart,yend,mxout,myout);
BeamProp=BeamProp;




