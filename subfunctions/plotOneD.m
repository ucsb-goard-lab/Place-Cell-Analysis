function [] = plotOneD(lap_by_lap_activity,valid_PCs,num_envs)
%PLOToneD Plots place cell responses by environment.
%
%   Input:
%   lap_by_lap_activity: a cell array where each cell contains a 3D array
%                        of lap by lap activity for each environment.
%   valid_PCs: a vector containing the indices of valid place cells.
%   num_envs: the number of environments.
%
%   This function plots the responses of place cells across different
%   environments. The responses are plotted 25 at a time, and the user can
%   press any key to see the next 25 cells. The responses are normalized
%   from 0 to 1 and the standard error across laps is also plotted.
%
%   The function first reshapes the lap by lap activity into a 4D array of
%   dimensions laps x bins x cells x envs, where the remainder is NaNs for
%   averaging. It then plots the responses of each cell across the
%   environments, with the responses in each environment being rescaled and
%   the standard error being plotted as well.
%
%   Note: This function requires the Statistics and Machine Learning
%   Toolbox to run the 'rescale' function.
%
%   Example:
%       plotOneD(lap_by_lap_activity,valid_PCs,3)
%
%   See also RESCALE, NANMEAN, NANSTD, SUBPLOT, PLOT.

% Plot PC responses 5 at a time by environment
figure;
sgtitle('Place Cell Responses by Environment')
color_opts = {'#7577CD','#FF6D02','#080957'}; % color options for plotting
color_opts_ste = [0.59,0.6,0.91;0.96,0.6,0.35;0.29,0.29,0.6]; % color options for plotting ste
disp(strcat('Total # place cells= ',num2str(length(valid_PCs))));
cell_count = 0; % start counting cells
n_bins = size(lap_by_lap_activity{1},2);

% reshape lap by lap activity to a 4D array of dims laps x bins x cells x
% envs (i.e. 15x72x2000x3) where remainder is NaNs (for averaging)
lap_nums = zeros(num_envs,1); % find # of laps run
for j =1:num_envs
    lap_nums(j,1) = size(lap_by_lap_activity{j},1);
end
laps = max(lap_nums);

lap_by_lap_array = NaN(laps,n_bins,size(lap_by_lap_activity{1},3),num_envs); 
for c = 1:size(lap_by_lap_activity{1},3) % for all cells
    for e = 1:num_envs
        curr_total_laps = NaN(laps,n_bins);
        curr_lap_activity = lap_by_lap_activity{e}(:,:,c);
        curr_total_laps(1:size(curr_lap_activity,1),:) = curr_lap_activity;
        lap_by_lap_array(:,:,c,e) = curr_total_laps;
    end
end

for ii = 1:ceil(length(valid_PCs)/25) % where 25 is # cells plotted per figure
    plot_count = 0; % subplot counter: should go from 1-25
    for i = 1:25
        plot_count = plot_count + 1; % augment subplot count
        cell_count = cell_count + 1; % augment cell count
%         smoothed_array = cat(3,activity_binned_1D_smoothed{1},activity_binned_1D_smoothed{2},activity_binned_1D_smoothed{3});
        if cell_count > length(valid_PCs)
            break % end plotting once last cell is reached
        end
        cell_resp = lap_by_lap_array(:,:,valid_PCs(cell_count),:);
        cell_resp = squeeze(cell_resp);
        avg_cell_resp = squeeze(mean(cell_resp,1,'omitnan'))'; % mean across laps
        ste_cell_resp = (squeeze(std(cell_resp,1,'omitnan'))/sqrt(laps))'; % get standard error across laps
%         cell_resp = smoothed_array(valid_PCs(cell_count), :, :);
%         cell_resp = squeeze(cell_resp)'; % array with DFF for one cell in each environment 
        [~,max_idx] = max(avg_cell_resp(1,:)); % location of max peak in env A
        %             centered_resp = zeros(num_envs,size(activity_binned_1D_smoothed{1,:}, 2)); % envs x bins
        % [~,max_idx] = max(cell_resp(:,:,1)); % location of max peak in env A
        %             dist2shift = 36-max_idx; % distance you'd have to circularly shift the vector to center it
        %             for jj = 1:num_envs
        %                 single_cent_resp = circshift(cell_resp(:,:,jj),dist2shift); % shift 1D response the requisite distance
        %                 centered_resp(jj,:) = single_cent_resp; % put shifted vector into 2D array
        %             end
        for j = 1:num_envs
            subplot(5,5,plot_count) % total = 25 cells
            x = 1:n_bins;
            curr_resp = rescale(avg_cell_resp(j,:)); % cell response for current environment normalized from 0 to 1
            rescale_factor = 0.5/max(avg_cell_resp(j,:));
            curve1 = curr_resp-ste_cell_resp(j,:)*rescale_factor; % lower bound ste
            curve2 = curr_resp+ste_cell_resp(j,:)*rescale_factor; % upper bound ste
            a = plot(x, curve1, 'r', 'LineWidth', 2);
            a.Color = color_opts_ste(j,:);
            hold on;
            b = plot(x, curve2, 'b', 'LineWidth', 2);
            b.Color = color_opts_ste(j,:);
            x2 = [x, fliplr(x)];
            inBetween = [curve1, fliplr(curve2)];
            fill(x2, inBetween, color_opts_ste(j,:));
            hold on
            p = plot(x, curr_resp, 'LineWidth', 1.5);
            p.Color = color_opts{j};
            if j == 1
                xline(max_idx) % plot a vertical line in the middle of the plot, for centered: line at 36
                title(strcat('Cell #',num2str(valid_PCs(cell_count)))); % plot cell # as subtitle
                xlabel('bins')
                ylabel('DFF')
            end
            hold on
        end
        hold off
    end
    pause % hit any key to see next 25 cells
end