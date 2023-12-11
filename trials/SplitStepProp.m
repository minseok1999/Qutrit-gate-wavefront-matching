function BeamProp = SplitStepProp(Beam,KZ,distance)
%%  Propagate beam via FFT
% - "Beam" is the transverse structure
% - "KZ" is a grid of wave vector components in the propagation direction
%   (this is a k-space grid and is dependant on the wavelength)
%   (Note: This grid needs to be in the fftshifted form where the k-space 
%   origin is in the corners of the grid)
% - "distance" is the propagation distance we want to simulate the propagation over
% The output is the "Beam"'s transverse structure after a propagation "distance"

%{
lamda=0.808e-6; % Wavelength (m)
    k=2*pi/lamda;

    
    
% Other parameters
    % Simulation window size (bigger helps with unwanted "computational box reflections"
    % in split step propagation):
    WidthX = 5e-3; % SLM window size or arbitrary size (m)
    WidthY = 5e-3; % SLM window size or arbitrary size (m)
     PixelSize = 0.008e-3; % width = height(m)
    % Grid:
    nx = round(WidthX/PixelSize); % Amount of pixels (x-dir)
    ny = round(WidthY/PixelSize); % Amount of pixels (y-dir)
    x = (-nx/2:nx/2-1)/nx*WidthX;            
    y = (-ny/2:ny/2-1)/ny*WidthY; 
    [X,Y] = meshgrid(x,y); % Cartesian grid
    
   

    
F0=exp(j*k*distance).*exp(j.*k.*(X^2+Y^2)./2./distance)./(j*lamda.*distance);
BeamProp=F0.*fft2(Beam.*exp(j.*pi./lamda./distance.*(X.^2+Y.^2)));

%}


global w_in
global WidthX
global WidthY
global PixelSize
% Grid:
    nx = round(WidthX/PixelSize); % Amount of pixels (x-dir)
    ny = round(WidthY/PixelSize); % Amount of pixels (y-dir)
    x = (-nx/2:nx/2-1)/nx*WidthX;            
    y = (-ny/2:ny/2-1)/ny*WidthY; 
    [X,Y] = meshgrid(x,y); % Cartesian grid
    clear x y
    Aperture=sign(1-sign(X.^2+Y.^2-((2.5).*w_in).^2));
    Beam=Beam.*Aperture;


BeamFFT = fft2(Beam);
BeamK = BeamFFT.*exp(1i*(KZ*distance));    
BeamProp = ifft2(BeamK);

