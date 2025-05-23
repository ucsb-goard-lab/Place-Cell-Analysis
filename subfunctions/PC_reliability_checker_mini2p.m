function [ PCs, criteria ] = PC_reliability_checker_mini2p(data, bin_size, coords)
%% PC_reliability_checker
%-------------------------------------------------------------------------%
%   This function finds place cells using the following three criteria, as
%   defined by Zong et al., 2022:
%
%   (1) Spatial information score > the 95th percentile of a shuffled
%   distribution
%   (2) Pearson's correlation between the 1st and 2nd half of the session >
%   the 80th percentile of a shuffled distribution
%   (3) At least one place fied, defined as activity elevated 40% above the 
%   mean in at least 2 square bins
%
%   Inputs:
%   'activityMap' is a 3D array where each slice along the third dimension 
%    contains the binned activity data.
%   'save_flag' is a boolean flag indicating whether to save data to the 
%    current directory. If 'save_flag' is not provided or empty, it defaults to 1 (true).
%   'shuffle_data' is a boolean flag indicating whether to shuffle actual 
%    data. If 'shuffle_data' is not provided or empty, it defaults to 0 (false).
%
%   Outputs:
%   'PCs' is a vector of indices of the place cells.
%   'FWHM' is a vector of the full width at half maximum (FWHM) of the 
%    place cells.
%   'params' is a structure containing the parameters used in the function.
%
%   Written by MJG 06/18/2021 // Last updated by NSW 06/01/2023
%
%   See also UIGETFILE, IMPORTDATA, RANDPERM, CORRCOEF, RANKSUM, FIT, MEAN, STD, SAVEAS.
%-------------------------------------------------------------------------%
% Parameters
nShuffles = 1000;
activity = data.DFF_transients;
x = coords(:,1);
y = coords(:,2);
[nCells, nFrames] = size(activity);

% Binning
xEdges = floor(min(x)):bin_size:ceil(max(x));
yEdges = floor(min(y)):bin_size:ceil(max(y));
nBinsX = length(xEdges)-1;
nBinsY = length(yEdges)-1;
nBins = nBinsX * nBinsY;
[~, ~, xBin] = histcounts(x, xEdges);
[~, ~, yBin] = histcounts(y, yEdges);

% Occupancy
occupancyMap = accumarray([yBin(xBin>0 & yBin>0), xBin(xBin>0 & yBin>0)], 1, [nBinsY, nBinsX]);
occupancyFlat = reshape(occupancyMap, [nBins, 1]);
p_i = occupancyFlat / sum(occupancyFlat);

% Bin index for each frame
disp('Getting spatial information score...')
binIdx = zeros(nFrames, 1);
for i = 1:nFrames
    if xBin(i) > 0 && yBin(i) > 0
        binIdx(i) = sub2ind([nBinsY, nBinsX], yBin(i), xBin(i));
    end
end

% Spatial information
activityFlat = zeros(nBins, nCells);
for c = 1:nCells
    tempMap = accumarray(binIdx(binIdx>0), activity(c, binIdx>0), [nBins, 1]);
    r_i = zeros(nBins,1);
    valid = occupancyFlat > 0;
    r_i(valid) = tempMap(valid) ./ occupancyFlat(valid);
    activityFlat(:,c) = r_i;
end

spatialInfo = zeros(1, nCells);
for c = 1:nCells
    r_i = activityFlat(:,c);
    r_i(r_i==0) = eps;
    r_bar = sum(p_i .* r_i);
    spatialInfo(c) = sum(p_i .* (r_i / r_bar) .* log2(r_i / r_bar));
end

% Shuffled SI
shuffledSI = zeros(nShuffles, nCells);
for s = 1:nShuffles
    for c = 1:nCells
        shuffled_trace = circshift(activity(c,:), randi(nFrames));
        tempMap = accumarray(binIdx(binIdx>0), shuffled_trace(binIdx>0), [nBins, 1]);
        r_i = zeros(nBins, 1);
        valid = occupancyFlat > 0;
        r_i(valid) = tempMap(valid) ./ occupancyFlat(valid);
        r_i(r_i==0) = eps;
        r_bar = sum(p_i .* r_i);
        shuffledSI(s,c) = sum(p_i .* (r_i / r_bar) .* log2(r_i / r_bar));
    end
end

sigSI = spatialInfo > prctile(shuffledSI, 95);

% Split-half reliability
disp('Computing consistency...')
half = floor(nFrames / 2);
reliability = zeros(1, nCells);
shuffledCorr = zeros(nShuffles, nCells);

for c = 1:nCells
    % First half
    idx1 = binIdx(1:half);
    vals1 = activity(c, 1:half);
    valid1 = idx1 > 0;
    act1 = accumarray(idx1(valid1), vals1(valid1), [nBins,1]);

    % Second half
    idx2 = binIdx(half+1:end);
    vals2 = activity(c, half+1:end);
    valid2 = idx2 > 0;
    act2 = accumarray(idx2(valid2), vals2(valid2), [nBins,1]);

    % Normalize by occupancy
    r1 = zeros(nBins,1); r2 = zeros(nBins,1);
    r1(valid) = act1(valid);
    r2(valid) = act2(valid);
    reliability(c) = corr(r1, r2, 'rows', 'complete');

    % Shuffle
    for s = 1:nShuffles
        s1 = circshift(vals1, randi([10,half], 1));
        s2 = circshift(vals2, randi([10,half], 1));

        act1s = accumarray(idx1(valid1), s1(valid1), [nBins,1]);
        act2s = accumarray(idx2(valid2), s2(valid2), [nBins,1]);

        r1s = zeros(nBins,1); r2s = zeros(nBins,1);
        r1s(valid) = act1s(valid);
        r2s(valid) = act2s(valid);
        shuffledCorr(s, c) = corr(r1s, r2s, 'rows', 'complete');
    end
end

sigRel = reliability > prctile(shuffledCorr, 70);

% Place field contiguity
disp('Checking for place fields...')
hasPlaceField = false(1, nCells);
thresholdFrac = 0.4;

for c = 1:nCells
    rMap = reshape(activityFlat(:,c), [nBinsY, nBinsX]);
    rThresh = mean(rMap(:)) + thresholdFrac * mean(rMap(:));
    above = rMap > rThresh;

    % Check for 2x2 block
    fieldDetected = false;
    for i = 1:(nBinsY - 1)
        for j = 1:(nBinsX - 1)
            block = above(i:i+1, j:j+1);
            if all(block(:))
                fieldDetected = true;
                break
            end
        end
        if fieldDetected, break; end
    end
    hasPlaceField(c) = fieldDetected;
end

% Final criteria
criteria = false(nCells, 3);  % [spatialInfo, reliability, placeField]
criteria(:,1) = sigSI(:);
criteria(:,2) = sigRel(:);
criteria(:,3) = hasPlaceField(:);

PCs = sigSI & sigRel & hasPlaceField;
fprintf('Identified %d place cells out of %d.\n', sum(PCs), nCells);
end