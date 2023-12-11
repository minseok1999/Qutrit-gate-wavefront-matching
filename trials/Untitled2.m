matfile = 'unitary_trial';	
Hologram=load(fullfile('C:\Users\rmins\Desktop\MPLC\MPLC_system\WaveFrontMatching_Code\trials',matfile));	
Hologram=Hologram.Hologram;

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