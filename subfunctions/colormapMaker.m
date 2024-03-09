function cmap = colormapMaker(colors, scaling)
%-------------------------------------------------------------------------%
 %COLORMAPMAKER Creates a colormap from a given set of colors.
 %   cmap = COLORMAPMAKER(colors) creates a colormap that transitions
 %   between the colors specified in the Nx3 array 'colors'. Each row of
 %   'colors' specifies a color in RGB format. The colormap transitions
 %   from each color to the next in the order provided. The scaling
 %   between colors is linear by default.
 %
 %   cmap = COLORMAPMAKER(colors, scaling) allows for specification of
 %   the scaling between colors. 'scaling' can be either 'linear' or
 %   'log', specifying a linear or logarithmic scaling, respectively.
 %
 %   Example:
 %       cmap = colormapMaker([255,255,255;0,255,0]);
 %   This example creates a colormap that transitions from white to green.
 %
 %   Note:
 %   - The input colors can be in the range 0-1 or 0-255. If any color
 %     component is greater than 1, it is assumed that the colors are in
 %     the range 0-255, and they are rescaled to 0-1.
 %   - The output colormap 'cmap' is an Mx3 array suitable for use with
 %     the colormap function in MATLAB. M is determined based on the
 %     number of input colors.
 %
 %   See also COLORMAP, LINESPACE, LOGSPACE.
%-------------------------------------------------------------------------%
    % takes colors (given as an N x 3 array) and creates a colormap that goes from each color in the order provided
    % ex. for white to green, colors = [255,255,255;0,255,0];
    if nargin < 2 || isempty(scaling)
        scaling = 'linear';
    end

    if size(colors, 2) ~= 3
        error('Make sure your input colors are N x 3, where N is the number of colors');
    end

    if any(colors > 1)
        colors = colors/255; % supplied 0-255 rather than 0-1, so we rescale
    end 

    n_colors = floor(255/(size(colors, 1)-1));

    colormap_matrix = zeros(3, n_colors);
    for n = 1:size(colors, 1) - 1
        switch scaling
        case 'linear'
            colormap_matrix(:, (n - 1) * n_colors + 1 : n * n_colors) = cat(1, linspace(colors(n, 1), colors(n+1, 1), n_colors), linspace(colors(n, 2), colors(n+1, 2), n_colors), linspace(colors(n, 3), colors(n+1, 3), n_colors));
        case 'log'
            colormap_matrix(:,  (n - 1) * n_colors + 1 : n * n_colors) = cat(1, logspace(colors(n, 1), colors(n+1, 1), n_colors), logspace(colors(n, 2), colors(n+1, 2), n_colors), logspace(colors(n, 3), colors(n+1, 3), n_colors));
        end
    end
    cmap = rescale(colormap_matrix)';
end