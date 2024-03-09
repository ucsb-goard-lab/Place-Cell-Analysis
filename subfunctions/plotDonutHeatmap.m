function [] = plotDonutHeatmap(valid_PCs,activity_1D)
%-------------------------------------------------------------------------%
%   PLOTDONUTHEATMAP(valid_PCs, activity_1D) plots the place cell responses by angular location
%   on a donut-shaped heatmap for each place cell in three different environments.
%
%   Inputs:
%   'valid_PCs' is a vector of indices of the place cells.
%   'activity_1D' is a cell array where each cell contains a 2D array of the binned activity data for a place cell in a specific environment.
%
%   This function does not return any outputs. Instead, it generates a series of plots, one for each place cell. Each plot shows the place cell responses by angular location on a donut-shaped heatmap.
%
%   See also ADDPATH, FIGURE, SGTITLE, CAT, LENGTH, MESHGRID, ONES, FIND, POLY2MASK, SUBPLOT, IMAGESC, AXIS, COLORMAP, TITLE, STRCAT, PAUSE.
%-------------------------------------------------------------------------%
%% Plot Neurotar data as a function of location along the circular track
figure;
sgtitle('Place Cell Responses by Angular Location')
color_opts = {[117,119,205],[255,109,2],[8,9,87]};
cell_count = 0;
center = 180;
binsize = 5;
addpath('E:\Code\Plotting') % for colormapMaker

pc_activity = cat(3,activity_1D{1}(valid_PCs,:),activity_1D{2}(valid_PCs,:),...
    activity_1D{3}(valid_PCs,:)); % put place cell data into 3D array
n_env = length(activity_1D);

x = -1*center:center; % i.e. -180:180
y = -1*center:center;
[xx yy] = meshgrid(x,y);
donut = ones(size(xx));
donut((xx.^2+yy.^2)<center^2)=0;   % big circle: radius 180, center at the origin
donut((xx.^2+yy.^2)<80^2)=1;   % small circle: radius 80, center at the origin
for ii = 1:ceil(size(pc_activity,1)/5) % for all cells
    plot_count = 0;
    for i = 1:5 % for all 5 subplots
        cell_count = cell_count + 1; % from 1 to # place cells
        if cell_count > length(valid_PCs)
            break  % end plotting once last cell is reached
        end
        curr_cell = pc_activity(cell_count,:,:);
        for ee = 1:n_env
            plot_count = plot_count + 1; % from 1 to 15
            curr_cell_env = rescale(curr_cell(:,:,ee)); % normalized single env DFF
            triangle = zeros(361); % create blank triangle matrix
            for bb = 1:size(curr_cell,2) % for each bin
                baseang2 = (binsize*2)*bb; % angle from base vector at (180,180) to (360,180)
                baseang1 = baseang2-binsize*2; % angle from base vector at (180,180) to (360,180)
                theta1 = (center-baseang1)/2; % botton two angles of isoceles
                theta2 = (center-baseang2)/2;
                x2 = center+center*cosd(theta1); % second x vertex
                y2 = center+center*sind(theta1);
                x3 = center+center*cosd(theta2); % third x vertex
                y3 = center+center*sind(theta2);
                x = [center,x2,x3]; % x vertices of triangle
                y = [center,y2,y3]; % y vertices of triangle
                mask = poly2mask(x,y,size(triangle,1),size(triangle,2));
                triangle(mask) = curr_cell_env(bb); % fill area with specified value
            end

            triangle(find(donut)) = 0; % erase all areas inside small circle and outside big circle
            p = subplot(5,3,plot_count); % 5 cells x 3 envs
            imagesc(triangle)
            axis square
            cmap = colormapMaker([255,255,255;color_opts{1,ee}]); % make colormap from white to env color
            colormap(p,cmap)
            if ee == 2 % just for second enviornment
                title(strcat('Cell #',num2str(valid_PCs(cell_count)))); % plot cell # as subtitle
            end
        end
    end
    pause % press any key to advance
end
