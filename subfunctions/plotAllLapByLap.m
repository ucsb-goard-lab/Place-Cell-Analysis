function [] = plotAllLapByLap(lap_by_lap_activity,valid_PCs)
%-------------------------------------------------------------------------%
%   PLOTALLLAPBYLAP(lap_by_lap_activity, valid_PCs) plots the lap-by-lap activity
%   for each place cell in three different environments: A, B, and A'.
%
%   Inputs:
%   'lap_by_lap_activity' is a 3D array where each slice along the third dimension contains the binned activity data for a lap.
%   'valid_PCs' is a vector of indices of the place cells.
%
%   This function does not return any outputs. Instead, it generates a series of plots, one for each place cell. Each plot shows the lap-by-lap activity of the place cell in the three environments.
%
%   See also ADDPATH, LENGTH, FIGURE, SUBPLOT, IMAGESC, COLORMAP, TITLE, SGTITLE, MIN, MAX, COLORBAR, SET, CLIM, PAUSE.
%-------------------------------------------------------------------------%

addpath('E:\Code\Plotting')
n_cells = length(valid_PCs);
envA = lap_by_lap_activity{1,:}(:,:,valid_PCs); % just place cells from env
envB = lap_by_lap_activity{2,:}(:,:,valid_PCs);
envA2 = lap_by_lap_activity{3,:}(:,:,valid_PCs);
figure
for ii = 1:n_cells
    subplot(3,1,1)
    imagesc(envA(:,:,ii))
    c = colormapMaker([255,255,255;23,107,135]); % white to blue
    colormap(c)
    title('Environment A')
    subplot(3,1,2)
    imagesc(envB(:,:,ii))
    colormap(c)
    title('Environment B')
    subplot(3,1,3)
    imagesc(envA2(:,:,ii));
    colormap(c)
    title('Environment A''')
    sgtitle(strcat('Lap by Lap Activity: Cell #', num2str(valid_PCs(ii))))

    minColorLimit = min([min(min(envA(:,:,ii))),min(min(envB(:,:,ii))),min(min(envA2(:,:,ii)))]);
    maxColorLimit = max([max(max(envA(:,:,ii))),max(max(envB(:,:,ii))),max(max(envA2(:,:,ii)))]);
    c = colorbar;
    set(c,'Position',[0.93 0.168 0.022 0.7])  % attach colorbar to axis
    clim([minColorLimit,maxColorLimit]);        % set colorbar limits
    pause
end