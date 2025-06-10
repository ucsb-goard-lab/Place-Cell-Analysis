function out = nanGaussFilt2D(data, sigma)
    % Create Gaussian kernel
    kSize = 2 * ceil(3 * sigma) + 1; % usually ±3 sigma
    g = fspecial('gaussian', kSize, sigma);

    % Replace NaNs with 0 temporarily
    dataZeroed = data;
    dataZeroed(isnan(data)) = 0;

    % Convolve image and mask
    weights = double(~isnan(data));
    filteredData = imfilter(dataZeroed, g, 'replicate');
    normFactor = imfilter(weights, g, 'replicate');

    % Avoid division by zero
    normFactor(normFactor == 0) = NaN;
    out = filteredData ./ normFactor;
end