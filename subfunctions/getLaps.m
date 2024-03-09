function [lap_start_times, lap_end_times, n_laps] = getLaps(binned_angles,...
    trialID, min_lap_length, max_avg_lap_speed, T, n_bins, floating, plot_flag)
%GETLAPS Finds the start and end times of laps in a mouse's path based on binned angles.
%   INPUT:
%   'binned_angles' is a vector of the mouse's angular position at each time point.
%   'trialID' is a string identifier for the trial.
%   'min_lap_length' is the minimum number of time points for a lap.
%   'max_avg_lap_speed' is the maximum average speed for a lap.
%   'T' is the total number of time points.
%   'n_bins' is the number of bins used for binning the angles.
%   'floating' is a structure with a field 'speed', which is a vector of the mouse's speed at each time point.
%   'plot_flag' is a boolean flag indicating whether to plot the binned angles and lap times.
%
%   OUTPUT:
%   The function returns 'lap_start_times' and 'lap_end_times', which are vectors of the start and end times of each lap,
%   and 'n_laps', which is the total number of laps.
%
%   This function follows a two-step processing. First, it finds all the times where
%   the mouse passes through the position it started the session at after
%   not being at that position for a while. Second, it makes sure that the
%   average speed during that lap was not too fast.
%
%   UPDATED 08/15/23 NSW

if strcmp(trialID, 'TRIAL_5') || strcmp(trialID, 'TRIAL_7') || isempty(trialID) % for trials 5 and 7 recordings: different neurotar version
    lap_start_times = find((binned_angles(2:end) - binned_angles(1:(end - 1)) < -30)) + 1;
    backward_lap_times = find((binned_angles(2:end) - binned_angles(1:(end - 1)) > 40)); % find when mouse backtracks across lap line
else % for trial 6 recordings: different neurotar version
    lap_start_times = find((binned_angles(2:end) - binned_angles(1:(end - 1)) > 30)) + 1;
    backward_lap_times = find((binned_angles(2:end) - binned_angles(1:(end - 1)) < -50));
end
lap_end_times = [lap_start_times - 1, T];
lap_start_times = [1,lap_start_times]; % add back in first position to start times

% remove times when mouse ran backwards
for bb = 1:length(backward_lap_times)
    curr_back_time = backward_lap_times(bb);
    % if there's not a whole lap in between backward point and last lap start time, delete lap start time
    early_laps = lap_start_times(lap_start_times < curr_back_time);
    last_start = max(early_laps);
    last_lap_angles = binned_angles(last_start+1:curr_back_time-2); % 2 & 1 for buffer
    curr_start = find(lap_start_times==last_start);
    if ~(max(last_lap_angles) - min(last_lap_angles) > 50) % if there's not a lap in between
        lap_start_times(curr_start) = []; % delete backwards lap start time
    elseif ismember(curr_back_time,lap_start_times) % if it got included as a regular start time
        lap_start_times(find(lap_start_times==curr_back_time)) = [];
    elseif isempty(last_lap_angles)
        lap_start_times(curr_start) = []; % if start time is within two positions of backward time
    end
    % if there's a plateau between backward point and next lap end time, delete lap end time
    late_laps = lap_end_times(lap_end_times > curr_back_time);
    next_end = min(late_laps);
    next_lap_angles = binned_angles(curr_back_time+2:next_end-1); % 2 & 1 for buffer
    if max(next_lap_angles) - min(next_lap_angles) < 10 % if there's a plateau in between
        lap_end_times(find(lap_end_times==next_end)) = []; % delete backwards lap start time
    elseif isempty(next_lap_angles) % if end time is within two positions of backward time
        lap_end_times(find(lap_end_times==next_end)) = [];
    else % if the lap end isn't deleted, add start time back in (very rare)
        if curr_start == 1
            lap_start_times = [last_start,lap_start_times];
        elseif curr_start > length(lap_start_times)
            lap_start_times = [lap_start_times,last_start];
        else
            lap_start_times = [lap_start_times(1:curr_start-1),last_start,...
                lap_start_times(curr_start:end)];
        end
    end
end

% remove any double end or start times
% todelete_e = [];
% for ee = 1:length(lap_end_times)-1
%     e_diff = lap_end_times(ee+1) - lap_end_times(ee);
%     if e_diff < 2
%         todelete_e = ee;
%     end
% end
% lap_end_times(todelete_e) = [];
% todelete_s = [];
% for ss = 1:length(lap_start_times)-1
%     s_diff = lap_start_times(ss+1) - lap_start_times(ss);
%     if s_diff < 2
%         todelete_s = ss;
%     end
% end
% lap_start_times(todelete_s) = [];


% remove laps that do not satisfy parameters
to_delete = [];
for ll = 1:length(lap_start_times)
    lap_length = lap_end_times(ll) - lap_start_times(ll);

    if lap_length < min_lap_length || abs(max(binned_angles(lap_start_times(ll):lap_end_times(ll)) - min(binned_angles(lap_start_times(ll):lap_end_times(ll))))) < (n_bins / 2)
        to_delete = [to_delete, ll];
    end
end
lap_start_times(to_delete) = [];
lap_end_times(to_delete) = [];

avg_speed = zeros(1, length(lap_start_times));
for ll = 1:length(lap_start_times)
    avg_speed(ll)= mean(floating.speed(lap_start_times(ll):lap_end_times(ll)));
end

lap_start_times(avg_speed > max_avg_lap_speed) = [];
lap_end_times(avg_speed > max_avg_lap_speed) = [];

n_laps = length(lap_start_times);

% plot start and end points on laps
if plot_flag
    envs = {'Env A','Env B','Env A'''};
    figure
    for ii = 1:length(data.binned_angles)
        binned_angles = data.binned_angles{ii};
        lap_start_times = data.lap_times{ii}(:,1);
        lap_end_times = data.lap_times{ii}(:,2);
        start_end_times = NaN(size(binned_angles));
        if strcmp(trialID,'TRIAL_6')
            start_end_times(lap_start_times) = 72;
            start_end_times(lap_end_times) = 1;
        else
            start_end_times(lap_start_times) = 1;
            start_end_times(lap_end_times) = 72;
        end
        x = 1:length(binned_angles);
        % figure
        subplot(3,1,ii)
        plot(binned_angles,'LineWidth',2)
        hold on
        scatter(x,start_end_times,'MarkerEdgeColor','m')
        title(envs{ii},'FontSize',12)
        xlabel('Time')
        ylabel('bin')
    end
end