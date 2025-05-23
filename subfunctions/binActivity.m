function [ activityMap ] = ...
    binActivity(data, coords, bin_size)
%-------------------------------------------------------------------------%
%   binActivity Bins the arena data by pixels.
%
%   Inputs:
%   'data' is a structure containing the suite2p activity data, with fields 
%   'DFF' and 'DFF_transients'.
%   'coords' is a two-column matrix containing the x and y coordinates of
%    the mouse.
%   'bin_size' is the N x N size of the bins.
%   'framesize' is the 3 dimensional size of the camera frames.
%
%   Outputs:
%   'activityMap' is a 3D array of binned activity.
%
%   Written by WTR 09/16/2020 // Last updated by NSW 05/15/2025
%-------------------------------------------------------------------------%
disp('Binning cell activity...')

%% Globals
N = size(data.DFF, 1);
activity = data.DFF_transients;

% Extract individual coordinates
x = coords(:,1);
y = coords(:,2);

% Define bin edges based on the full arena extent
xEdges = floor(min(x)):bin_size:ceil(max(x));
yEdges = floor(min(y)):bin_size:ceil(max(y));
numBinsX = length(xEdges) - 1;
numBinsY = length(yEdges) - 1;

% Compute 2D histogram of occupancy
[occupancyCounts, ~, xBin] = histcounts(x, xEdges);
[~, ~, yBin] = histcounts(y, yEdges);

% Initialize
activityMap = zeros(numBinsY, numBinsX, N);
occupancyMap = zeros(numBinsY, numBinsX);

% Compute occupancy map
for i = 1:length(x)
    xi = xBin(i);
    yi = yBin(i);
    if xi > 0 && yi > 0
        occupancyMap(yi, xi) = occupancyMap(yi, xi) + 1;
    end
end

% Compute activity maps
for cellIdx = 1:N
    for i = 1:length(x)
        xi = xBin(i);
        yi = yBin(i);
        if xi > 0 && yi > 0
            activityMap(yi, xi, cellIdx) = activityMap(yi, xi, cellIdx) + activity(cellIdx, i);
        end
    end
end

% Normalize by occupancy
for cellIdx = 1:N
    cellMap = activityMap(:,:,cellIdx);
    normMap = cellMap ./ occupancyMap;
    normMap(occupancyMap == 0) = NaN; % avoid dividing by zero
    activityMap(:,:,cellIdx) = normMap;
end


% nCells = size(activityMap, 3);  % number of cells
% 
% for cellIdx = 1:nCells
%     imagesc(activityMap(:,:,cellIdx));
%     axis equal tight;
%     colorbar;
%     title(sprintf('Cell %d Activity Map', cellIdx));
% 
%     % Pause or wait for user input
%     pause(0.5);
% end