%% This subscript adds a -pi to pi colorbar to the latest subplot
originalSize = get(gca, 'Position');
cbr = colorbar;
cbr.Ticks = [-3.1 0 3.1];
cbr.TickLabelInterpreter = 'latex';
cbr.TickLabels = {'$-\pi$', '$0$', '$\pi$'};
cbr.FontSize = 26;
cbr.FontName = 'Times New Roman';
set(gca, 'Position', originalSize);