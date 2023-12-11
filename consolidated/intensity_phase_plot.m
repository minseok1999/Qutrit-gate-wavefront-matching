function intensity_phase_plot(Mode_Efield)
%% Intensity_phase_plot
% Displays the transverse electric field's phase with a colormap and 
% the intensity overlayed as alpha (transparency) values, along with a 
% black background.

    Int = abs(Mode_Efield).^2;
    Int = Int/max(max(Int));
    Ang = angle(Mode_Efield);
    
    % Colormap for phase:
    colormap(hsv)
    imag = imagesc(Ang,[-pi pi]);
    
    %colormap(hsv) % Optional colormap 
    set(imag, 'AlphaData', Int) %overlay with intensity
    set(gca,'color',[0 0 0]) %set background black
    
    % Get rid of axes:
    set(gca,'xtick',[])
    set(gca,'ytick',[]) 

end