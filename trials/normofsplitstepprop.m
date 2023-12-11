
global lamda k
lamda=0.808e-6; % Wavelength (m)
    k=2*pi/lamda;
    
   
    PixelSize = 0.008e-3; % width = height(m)
    w_out = 0.9e-3; % Output beam radius (m)
    w_in = 0.9e-3; % Input beam radius (m)
    PropDist = 800e-3; % Propagation between phase screens (m)

% Other parameters
    % Simulation window size (bigger helps with unwanted "computational box reflections"
    % in split step propagation):
    WidthX = 5e-3; % SLM window size or arbitrary size (m)
    WidthY = 5e-3; % SLM window size or arbitrary size (m)

% Grid:
    nx = round(WidthX/PixelSize); % Amount of pixels (x-dir)
    ny = round(WidthY/PixelSize); % Amount of pixels (y-dir)
    x = (-nx/2:nx/2-1)/nx*WidthX;            
    y = (-ny/2:ny/2-1)/ny*WidthY; 
    [X,Y] = meshgrid(x,y); % Cartesian grid
    clear x y

    % Grid in cylindrical coordinates:
    Rad = sqrt(X.^2+Y.^2);
    Angle = angle(X+1i.*Y)+pi; % Matrix with all angles starting left-center

    
% K-space grid used in split step propagation:
    % Move the origin to match the SplitStepProp-function:
    kx = (mod(nx/2:nx+nx/2-1,nx)-(nx/2))*2*pi/WidthX;
    ky = (mod(ny/2:ny+ny/2-1,ny)-(ny/2))*2*pi/WidthY;
    [KX,KY] = meshgrid(kx,ky);
    % Matrix for k-space propagation direction components:
    KZ = sqrt((2*pi/lamda)^2-(KX.^2+KY.^2)); 
    clear kx ky KX KY
    
  gg = ((-1+i).*GenModesLG([1 0], w_out, Rad, Angle)+(1+i).*GenModesLG([-1 0], w_out, Rad, Angle)+(0).*GenModesLG([0 0], w_out, Rad, Angle))./2;
    
  norm1=(sum(abs(gg),'all')^2)
  
  ggg=SplitStepProp(gg,KZ,PropDist);
  
  norm2=(sum(abs(ggg),'all')^2)
  