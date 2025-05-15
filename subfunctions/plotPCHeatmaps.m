function [] = plotPCHeatmaps(activityMap, valid_PCs)

% Parameters
cellsPerFigure = 20;
cellIndices = find(valid_PCs);  % indices of identified place cells
nPlaceCells = length(cellIndices);
figure;
for batchStart = 1:cellsPerFigure:nPlaceCells
    batchEnd = min(batchStart + cellsPerFigure - 1, nPlaceCells);

    for i = batchStart:batchEnd
        subplot(4, 5, i - batchStart + 1); % 4 rows x 5 cols layout
        imagesc(activityMap(:,:,cellIndices(i)));
        axis equal tight off;
        colorbar;
        title(sprintf('Cell %d', cellIndices(i)));
    end

    sgtitle(sprintf('Place Cells %d–%d of %d', ...
        cellIndices(batchStart), cellIndices(batchEnd), nPlaceCells));

    pause;
end