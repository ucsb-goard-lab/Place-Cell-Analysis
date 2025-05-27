function [] = plotPCHeatmaps(activityMap, valid_PCs)
%-------------------------------------------------------------------------%
%   Plots smoothed heatmaps of place cell activity.
%   Applies Gaussian smoothing to each heatmap prior to display.
%
%   Inputs:
%       activityMap - 3D matrix (x, y, cell)
%       valid_PCs - logical or index array of place cell IDs
%
%   Written by NSW, updated 5/23/25 to include Gaussian smoothing
%-------------------------------------------------------------------------%

% Parameters
cellsPerFigure = 20;
sigma = 1;  % standard deviation of the Gaussian filter

cellIndices = find(valid_PCs);  % indices of identified place cells
nPlaceCells = length(cellIndices);

figure;
for batchStart = 1:cellsPerFigure:nPlaceCells
    batchEnd = min(batchStart + cellsPerFigure - 1, nPlaceCells);

    for i = batchStart:batchEnd
        subplot(4, 5, i - batchStart + 1); % 4 rows x 5 cols layout
        
        % Apply Gaussian smoothing to each cell's heatmap
        smoothedMap = imgaussfilt(activityMap(:,:,cellIndices(i)), sigma);

        imagesc(smoothedMap);
        axis equal tight off;
        colorbar;
        title(sprintf('Cell %d', cellIndices(i)));
    end

    sgtitle(sprintf('Place Cells %d–%d of %d', ...
        cellIndices(batchStart), cellIndices(batchEnd), nPlaceCells));

    pause;
end
