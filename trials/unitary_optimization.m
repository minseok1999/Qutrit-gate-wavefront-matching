%% Wavefront Matching script
% The script takes one mode as an input and one as an 
% output. After this, it propagates this input forwards, and the output 
% backwards through the system. Then, it calculates the overlaps of these 
% propagated modes at every phase screen (PhaseScreenNum) separated by a 
% set distance (PropDist). The structures of the phase screens are then 
% updated, so that the input modes eventually convert to the output modes,
% in the system, after a sufficient number of updates.

% The code was written for and utilized in the work described in:
% Hiekkam�ki, M., Prabhakar, S., & Fickler, R.
% "Near-perfect measuring of full-field transverse-spatial modes of light "
% Optics Express 27 (22), 31456 (2019) 

% The general idea is based on the article:
% Fontaine, N. K., Ryf, R., Chen, H., Neilson, D., & Carpenter, J.  
% "Design of high order mode-multiplexers using multiplane light conversion." 
% 2017 European Conference on Optical Communication (ECOC). IEEE, 2017.

% This version of the code requires:
%
% -> included in the zip-file
% AddPhaseColorbar.m (subscript)
% GenModesLG.m (function)
% intensity_phase_plot.m (function)
% SplitStepProp.m (function)
% update_PhaseScreens.m (subscript)
% 
% -> needs to be downloaded from Mathworks:
% LaguerreGen.m - function by Mattthias Trampisch (2019). 
%       Generalized Laguerre polynomial 
%       (https://www.mathworks.com/matlabcentral/fileexchange/15916-generalized-laguerre-polynomial) 
%       MATLAB Central File Exchange


% Created by: Markus Hiekkam�ki, Shashi Prabhakar, Robert Fickler, 
% markus.hiekkamaki@tuni.fi, shashi.sinha@tuni.fi, robert.fickler@tuni.fi, 
% Last Updated: 11.12.2019

clear all;clc;close all;format compact

%% Constants and Parameters:

% Regarding optimization:
    % Convergence accuracy (number of decimals):
    convLim = 4;
    % Number of previous iterations to consider when checking for convergence:
    convValue = 3;
    % Maximum number of iterations:
    IterMax = 1000;
    
% SLM specifications:
    PixelSize = 0.008e-3; % width = height(m)
    PhaseScreenNum = 3; % Number of Phase screens used in optimization
    MaxPhaseValue = 255; % Number of discrete phase shift values in the 2pi range

% Setup parameters:

global lamda k
lamda=0.808e-6; % Wavelength (m)
    k=2*pi/lamda;
    
    w_out = 0.9e-3; % Output beam radius (m)
    w_in = 0.9e-3; % Input beam radius (m)
    PropDist = 800e-3; % Propagation between phase screens (m)

% Other parameters
    % Simulation window size (bigger helps with unwanted "computational box reflections"
    % in split step propagation):
    WidthX = 5e-3; % SLM window size or arbitrary size (m)
    WidthY = 5e-3; % SLM window size or arbitrary size (m)
    
%% Flags controlling optimization:
   
% When are the phase screens updated: (Both can be used at the same time)
% (This can slightly affect which phase screen removes phase singularities
% [first or last usually])
    % Update when propagating backwards:
    BackUpdate = true;
    % Update when propagating forwards:
    ForwardUpdate = false;
    
% Display intermediate modes after optimization:
    DispInterim = true;
% Save results to the files specified below (false, if you don't want to save):
    SaveFlag = true;
    SaveBitmap = false;
    
%% Saving parameters:

% Name of folder for saving Bitmap holograms:
    picLoc = "Test_holos";
% Name of .mat file to save data to:
    FileName = "trial_superposition_10input.mat";

%% Input and output modes:

% Input: (Only LG modes here)
    % First colum is OAM (azimuthal) index, second is radial index
    % Example: 
    % InputModes=[-1 0]; Here the mode is an OAM=-1 mode 
    % without radial structure    
    InputModes1 = [1 0];
    InputModes2 = [-1 0];
% Ouput: 
    OutputModes1 = [1 0];
    OutputModes2 = [-1 0];

%% (NO FURTHER INPUT PARAMETERS FROM THIS POINT ONWARDS)     

%% Calculate rest of the parameters from inputs:

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
    
%% Calculate used modes:
fprintf('Calculating modes \n \n')

ModesIn1 = GenModesLG(InputModes1, w_in, Rad, Angle);
ModesOut1 = (GenModesLG([1 0], w_out, Rad, Angle)+GenModesLG([-1 0], w_out, Rad, Angle))./sqrt(2);
ModesIn2 = GenModesLG(InputModes2, w_in, Rad, Angle);
ModesOut2 = (GenModesLG([1 0], w_out, Rad, Angle)-GenModesLG([-1 0], w_out, Rad, Angle))./sqrt(2);

% Display modes:
InitialmodeFig = figure;
% Input modes1:               
subplot(2, 1, 1);
intensity_phase_plot(ModesIn1)
title('Input mode1')
AddPhaseColorbar

% Output modes1:     
subplot(2, 1, 2);
intensity_phase_plot(ModesOut1)
title('Output mode1')
AddPhaseColorbar

% Display modes:
InitialmodeFig = figure;

% Input modes:               
subplot(2, 1, 1);
intensity_phase_plot(ModesIn2)
title('Input mode2')
AddPhaseColorbar

% Output modes:     
subplot(2, 1, 2);
intensity_phase_plot(ModesOut2)
title('Output mode2')
AddPhaseColorbar

%% Propagate the input and output before starting optimization:
fprintf('Initial propagation \n \n')

% Initialize datastructures:
    Beam1 = zeros(nx,nx,PhaseScreenNum + 1); 
    % The (+ 1) Beam structure is used to store the simulated output that 
    % is compared with the wanted output.
    BeamBack1 = zeros(nx,nx,PhaseScreenNum);
    
    Hologram = zeros(nx,nx,PhaseScreenNum);
    
% Propagating the input forwards:
for PhScrInd = 1:1:PhaseScreenNum % For all phase screens
      
    if PhScrInd == 1 % = The Initial mode
        Beam1(:,:,1) = ModesIn1;
    end
   
    
    % Propagate beam from previous phase screen while applying the
    % phase screen's hologram
    % (hologram's structure doesn't matter if it's initially just zeros):
    Beam1(:,:,PhScrInd+1) = Beam1(:,:,PhScrInd).*exp(1i*(Hologram(:,:,PhScrInd)));
    Beam1(:,:,PhScrInd+1) = SplitStepProp(Beam1(:,:,PhScrInd+1),KZ,PropDist);
end

% Initialize datastructures:
    Beam2 = zeros(nx,nx,PhaseScreenNum + 1); 
    % The (+ 1) Beam structure is used to store the simulated output that 
    % is compared with the wanted output.
    BeamBack2 = zeros(nx,nx,PhaseScreenNum);
    
    for PhScrInd = 1:1:PhaseScreenNum % For all phase screens
      
    if PhScrInd == 1 % = The Initial mode
        Beam2(:,:,1) = ModesIn2;
    end
   
    
    % Propagate beam from previous phase screen while applying the
    % phase screen's hologram
    % (hologram's structure doesn't matter if it's initially just zeros):
    Beam2(:,:,PhScrInd+1) = Beam2(:,:,PhScrInd).*exp(1i*(Hologram(:,:,PhScrInd)));
    Beam2(:,:,PhScrInd+1) = SplitStepProp(Beam2(:,:,PhScrInd+1),KZ,PropDist);
end
% Then propagate the output backwards:
for PhScrInd=PhaseScreenNum:-1:2 % For all phase screens
    
    if PhScrInd == PhaseScreenNum % = For final output mode
        
        % We invert the phase to get backwards propagation when
        % using the same Split step propagation function.
        % Note: the wanted output is simulated to be "formed" one
        % "PropDist" away from the last hologram
        BeamBack1(:,:,PhaseScreenNum) = SplitStepProp(ModesOut1,KZ,PropDist);
    end
    
    % Propagate from previous phase screen while applying the phase
    % screen hologram:
    BeamBack1(:,:,PhScrInd-1) = BeamBack1(:,:,PhScrInd).*exp(1i*(Hologram(:,:,PhScrInd)));
    BeamBack1(:,:,PhScrInd-1) = SplitStepProp(BeamBack1(:,:,PhScrInd-1),KZ,PropDist);
end
    % Then propagate the output backwards:
for PhScrInd=PhaseScreenNum:-1:2 % For all phase screens
    
    if PhScrInd == PhaseScreenNum % = For final output mode
        
        % We invert the phase to get backwards propagation when
        % using the same Split step propagation function.
        % Note: the wanted output is simulated to be "formed" one
        % "PropDist" away from the last hologram
        BeamBack2(:,:,PhaseScreenNum) = SplitStepProp(ModesOut2,KZ,PropDist);
    end
    
    % Propagate from previous phase screen while applying the phase
    % screen hologram:
    BeamBack2(:,:,PhScrInd-1) = BeamBack2(:,:,PhScrInd).*exp(1i*(Hologram(:,:,PhScrInd)));
    BeamBack2(:,:,PhScrInd-1) = SplitStepProp(BeamBack2(:,:,PhScrInd-1),KZ,PropDist);
end
%% Creating a stop-button to manually stop the optimization at any point:
% Pressing this stop-button only stops the optimization and the results 
% will still be plotted. 
StopFig = figure;
set(StopFig,'Position',[300 300 100 50]); %figure appears at initial coordinates
set(StopFig,'menubar','none','units','pixels');
StopVar = uicontrol('Style', 'PushButton', 'String', 'Stop WFM', 'Callback', 'delete(gcbo)', 'Position', [1 1 100 50]);
drawnow

%% WaveFront Matching:
fprintf('WaveFront Matching \n \n')

% Initialize figures and parameters for observing the convergence:
    IterationCount = 0; % Tracks the iteration count of the process 
    % To store the visibility and overlap of optimized modes:
    resultoverlap1 = NaN(1, IterMax);
    
    % Initialize figure for tracking visibility and overlap:
    OptimFig = figure;
    OptimAx = subplot(1,1,1);
    line1 = plot(OptimAx,resultoverlap1);
    legend('Overlap', 'Location', 'southeast')
    ylabel('Overlap')
    xlabel('Iteration')
    grid on
    drawnow
    
      % To store the visibility and overlap of optimized modes:
    resultoverlap2 = NaN(1, IterMax);
    
    % Initialize figure for tracking visibility and overlap:
    OptimFig = figure;
    OptimAx = subplot(1,1,1);
    line2 = plot(OptimAx,resultoverlap2);
    legend('Overlap', 'Location', 'southeast')
    ylabel('Overlap')
    xlabel('Iteration')
    grid on
    drawnow
    
% Start optimizing:
% Until iterations reach IterNum or the stop button is pressed:
while IterationCount <= IterMax && ishandle(StopVar)
    IterationCount = IterationCount + 1;
    fprintf(strcat("Current iteration is ", num2str(IterationCount), "\n"))
    
% Propagate and update in the forwards direction
    for PhScrInd = 1:PhaseScreenNum
        
        if ForwardUpdate
            update_Phasescreens % Subscript to update current phase screen
        end
        
        % Modulate beam with new hologram:
        % Imprint hologram on forward propagating beam
        Beam1(:,:,PhScrInd+1) = Beam1(:,:,PhScrInd).*exp(1i*(Hologram(:,:,PhScrInd)));
        % Propagate to next hologram (or to the end)
        Beam1(:,:,PhScrInd+1) = SplitStepProp(Beam1(:,:,PhScrInd+1),KZ,PropDist);
            
    end
    
     for PhScrInd = 1:PhaseScreenNum
        
        if ForwardUpdate
            update_Phasescreens % Subscript to update current phase screen
        end
        
        % Modulate beam with new hologram:
        % Imprint hologram on forward propagating beam
        Beam2(:,:,PhScrInd+1) = Beam2(:,:,PhScrInd).*exp(1i*(Hologram(:,:,PhScrInd)));
        % Propagate to next hologram (or to the end)
        Beam2(:,:,PhScrInd+1) = SplitStepProp(Beam2(:,:,PhScrInd+1),KZ,PropDist);
            
     end
    
     
% Propagate and update in the backwards direction:
    for PhScrInd = PhaseScreenNum:-1:1
        
        if BackUpdate
            update_Phasescreens % Subscript to update current phase screen
        end
        
        % Modulate backwards beam with new hologram:
        if PhScrInd > 1
            BeamBack1(:,:,PhScrInd-1) = BeamBack1(:,:,PhScrInd).*exp(1i*Hologram(:,:,PhScrInd));
            % Propagate backwards (already conjugated in the first
            % propagations)
            BeamBack1(:,:,PhScrInd-1) = SplitStepProp(BeamBack1(:,:,PhScrInd-1),KZ,PropDist);
        end
        
          % Modulate backwards beam with new hologram:
        if PhScrInd > 1
            BeamBack2(:,:,PhScrInd-1) = BeamBack2(:,:,PhScrInd).*exp(1i*Hologram(:,:,PhScrInd));
            % Propagate backwards (already conjugated in the first
            % propagations)
            BeamBack2(:,:,PhScrInd-1) = SplitStepProp(BeamBack2(:,:,PhScrInd-1),KZ,PropDist);
        end
        
    end
    
    for PhScrInd = PhaseScreenNum:-1:1
        
        if BackUpdate
            update_Phasescreens % Subscript to update current phase screen
        end
        
    end
    
    
% Modulate forwards beam with new hologram:
    for PhScrInd=1:PhaseScreenNum
        % imprint hologram on forward beam
        Beam1(:,:,PhScrInd+1) = Beam1(:,:,PhScrInd).*exp(1i*(Hologram(:,:,PhScrInd)));
        % propagate to next hologram
        Beam1(:,:,PhScrInd+1) = SplitStepProp(Beam1(:,:,PhScrInd+1),KZ,PropDist);
    end
    for PhScrInd=1:PhaseScreenNum
        % imprint hologram on forward beam
        Beam2(:,:,PhScrInd+1) = Beam2(:,:,PhScrInd).*exp(1i*(Hologram(:,:,PhScrInd)));
        % propagate to next hologram
        Beam2(:,:,PhScrInd+1) = SplitStepProp(Beam2(:,:,PhScrInd+1),KZ,PropDist);
    end
% Check overlap and visibility, and plot current values
        
    % Check overlap between transformations and wanted outputs:
    Overlap_to_output1 = real(abs(sum(sum(conj(Beam1(:,:,PhaseScreenNum+1)).*ModesOut1))).^2);    
    resultoverlap1(IterationCount) = Overlap_to_output1;
    % Check overlap between transformations and wanted outputs:
    Overlap_to_output2 = real(abs(sum(sum(conj(Beam2(:,:,PhaseScreenNum+1)).*ModesOut2))).^2);    
    resultoverlap2(IterationCount) = Overlap_to_output2;
    
    try
        line1.YData = resultoverlap1;
        drawnow
    catch ME
        fprintf("Figure for optimization was closed.")
    end
         try
        line2.YData = resultoverlap2;
        drawnow
    catch ME
        fprintf("Figure for optimization was closed.")
         end
    
    % Display values related to converging:
    if IterationCount > convValue+1
        % Latest values:
        LastOverl = round(mean(resultoverlap1(IterationCount-convValue:IterationCount),'omitnan'), convLim);
        disp(strcat(['Current Overlap: ', num2str(resultoverlap1(IterationCount)), '; Mean of previous: ', num2str(LastOverl)]))
             
        LastOver2 = round(mean(resultoverlap2(IterationCount-convValue:IterationCount),'omitnan'), convLim);
        disp(strcat(['Current Overlap: ', num2str(resultoverlap2(IterationCount)), '; Mean of previous: ', num2str(LastOver2)]))
       
        % Stop optimization if the values converged:
        if round(resultoverlap1(IterationCount), convLim) <= LastOverl & round(resultoverlap2(IterationCount), convLim) <= LastOver2
            close(StopFig)
        end
    end
    
end

% To make sure the stop button and OptimFig are closed:
try
    close(StopFig)
    close(OptimFig)
catch ME
    try 
        close(OptimFig)
    catch ME
    end
end

 
%% Display the results

% Display the phase screens:
figure
for PhScrInd=1:PhaseScreenNum    
    subplot(1,PhaseScreenNum,PhScrInd);
    % Colormap:
    colormap(hsv)
    imagesc(Hologram(:,:,PhScrInd));
    % Get rid of axes:
    set(gca,'xtick',[])
    set(gca,'ytick',[])

    title(strcat(['Hologram', num2str(PhScrInd)]))   
end
% Colorbar:
originalSize = get(gca, 'Position');
cbr = colorbar;
cbr.Ticks = [0 3.14 6.25];
cbr.TickLabelInterpreter = 'latex';
cbr.TickLabels = {'$0$', '$\pi$', '$2\pi$'};
cbr.FontSize = 26;
cbr.FontName = 'Times New Roman';
set(gca, 'Position', originalSize);

% The beam structure in its intermediate stages:
if DispInterim
    figure
    for PhaseScrInd = 1:PhaseScreenNum + 1
        subplot(1, PhaseScreenNum + 1, PhaseScrInd)
        intensity_phase_plot(Beam1(:,:,PhaseScrInd))
        title(strcat(['Mode after hologram_1 ', num2str(PhaseScrInd-1)]))
    end
    AddPhaseColorbar
end


if DispInterim
    figure
    for PhaseScrInd = 1:PhaseScreenNum + 1
        subplot(1, PhaseScreenNum + 1, PhaseScrInd)
        intensity_phase_plot(Beam2(:,:,PhaseScrInd))
        title(strcat(['Mode after hologram_2 ', num2str(PhaseScrInd-1)]))
    end
    AddPhaseColorbar
end


% Plot the simulated conversion of input and output: 
figure
BeamBefore = ModesIn1;
subplot(2,1,1);
intensity_phase_plot(BeamBefore)
title('Initial mode')
AddPhaseColorbar
BeamAfter = Beam1(:,:,PhaseScreenNum+1);
subplot(2,1,2);
intensity_phase_plot(BeamAfter)
title('Simulated output')
AddPhaseColorbar

% Show final overlap:
fprintf(strcat("\n", "Final average overlap_1 is ", num2str(trace(Overlap_to_output1)), "\n"));


% Plot the simulated conversion of input and output: 
figure
BeamBefore = ModesIn2;
subplot(2,1,1);
intensity_phase_plot(BeamBefore)
title('Initial mode')
AddPhaseColorbar
BeamAfter = Beam2(:,:,PhaseScreenNum+1);
subplot(2,1,2);
intensity_phase_plot(BeamAfter)
title('Simulated output')
AddPhaseColorbar

% Show final overlap:
fprintf(strcat("\n", "Final average overlap_2 is ", num2str(trace(Overlap_to_output2)), "\n"));
    
%% Save data

if SaveFlag
   save(FileName,'Hologram') 
end
% Save bitmap holograms:
if SaveBitmap
    % Check if given folder exists
    if ~exist(picLoc, 'dir')
           mkdir(picLoc)
    end
    % Save holograms to the given folder with names "H1.bmp",...
    for i = 1:size(Hologram, 3)
        imwrite(mat2gray(Hologram(:,:,i),[0, 2*pi]), strcat(picLoc, 'H',num2str(i),'.bmp'))
    end
end