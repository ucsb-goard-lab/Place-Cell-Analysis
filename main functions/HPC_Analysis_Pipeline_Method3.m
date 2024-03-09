function [data,active_cells,spatial_info,activity_binned_1D_smoothed,...
    activity_binned_1D_sde_smoothed,lap_by_lap_activity,params,valid_PCs,...
    valid_SCs,sorted_PCs,PC_activity] = HPC_Analysis_Pipeline_Method3(num_envs,...
    data, floating, plot_donut, plot_laps, plot_oneD, plot_heatMap, save_flag)
%-------------------------------------------------------------------------%
%HPC_ANALYSIS_PIPELINE_METHOD3 Analyzes hippocampal place cell data.
%   This function is an updated version of the HPC analysis pipeline that raises the 
%   consistency threshold and uses Cohen's d for place cell (PC) thresholding. It 
%   deconcatenates the floating data in the Neurotar Extractor and breaks it down 
%   into the length of each environment. It calculates where to split files based on 
%   the number of frames and saves the data into split arrays where N = num_envs.
%
%   Inputs:
%   - num_envs: The number of environments. If not provided or empty, the default value is 3.
%   - data: The processed data. If not provided or empty, the function will import data from 'data_2_1.mat'.
%   - floating: The floating data for the NewNeurotarExtractor class. If not provided or empty, the function will import data from 'floating_2_1.mat'.
%   - plot_donut: A flag to plot average activity on donut heatmap. If not provided or empty, the default value is 0 (no plot).
%   - plot_laps: A flag to plot lap by lap activity for every cell (not just PCs). If not provided or empty, the default value is 0 (no plot).
%   - plot_oneD: A flag to plot PC responses for each environment. If not provided or empty, the default value is 0 (no plot).
%   - plot_heatMap: A flag to plot heatmaps for each environment. If not provided or empty, the default value is 0 (no plot).
%   - save_flag: A flag to save data. If not provided or empty, the default value is 1 (save data).
%
%   Outputs:
%   - data: The deconcatenated data, split into the length of each environment.
%   - active_cells: A cell array containing the indices of active cells for each environment.
%   - spatial_info: A cell array containing the spatial information for each environment.
%   - activity_binned_1D_smoothed: A cell array containing the smoothed 1D binned activity for each environment.
%   - activity_binned_1D_sde_smoothed: A cell array containing the smoothed and standard error 1D binned activity for each environment.
%   - lap_by_lap_activity: A cell array containing the lap by lap activity for each environment.
%   - params: A structure containing the parameters used in the function.
%   - valid_PCs: A vector containing the indices of all cells with a place field in at least one environment.
%   - valid_SCs: A vector containing the indices of all speed cells in at least one environment.
%   - sorted_PCs: A vector containing the indices of sorted place cells (only if plot_heatMap is true).
%   - PC_activity: A matrix containing the activity of place cells (only if plot_heatMap is true).
%
%   Written by NSW 05/24/2022 // Last updated by NSW 08/31/2023
%-------------------------------------------------------------------------%
if nargin < 1 || isempty(num_envs)
    num_envs = 3;
end
if nargin < 2 || isempty(data)
    data = importdata('data_1.mat');
end
if nargin < 3 || isempty(floating)
    floating = importdata('floating_1.mat'); % for NewNeurotarExtractor class
end
if nargin < 4 || isempty(plot_donut)
    plot_donut = 0; % plots average activity on donut heatmap
end
if nargin < 5 || isempty(plot_laps)
    plot_laps = 0; % plot lap by lap activity for every cell (not just pc's)
end
if nargin < 6 || isempty(plot_oneD)
    plot_oneD = 0; % plot PC responses for each environment
end
if nargin < 7 || isempty(plot_heatMap)
    plot_heatMap = 0; % plot heatmaps for each environment
end
if nargin < 8 || isempty(save_flag)
    save_flag = 1; % save data
end

% preallocate for speed
[PC_array, SC_array, spatial_info, activity_binned_1D_smoothed,...
    consistency, speed_binned,active_cells,spikes_speed_binned,...
    sde_spikes_speed_binned, lap_by_lap_activity, FWHM,...
    angle_bin_ids,lap_times,DFF_trans] = deal(cell(num_envs,1));

for ii = 1:num_envs
    % Display current iteration
    formatSpec = 'Calculating responses for environment %1.0f...';
    fprintf(formatSpec,ii);

    if num_envs > 1 && ii > 1
        %% Data
        % Activity: 'data_N.mat'
        activity_fspec = 'data_%d.mat';
        activity_fname = sprintf(activity_fspec,ii);
        data = importdata(activity_fname);
    end

    if num_envs == 1
        fall_fname = 'Fall.mat';
        Fall_File = importdata(fall_fname);
    else
        % 'Fall_N.mat'
        fall_fspec = 'Fall_%d.mat';
        fall_fname = sprintf(fall_fspec,ii);
        Fall_File = importdata(fall_fname);
    end
    data.DFF(Fall_File.iscell(:, 1) == 0, :) = []; % delete DFF traces and spikes from non-cells
    data.spikes(Fall_File.iscell(:, 1) == 0, :) = [];
    data.DFF_raw(Fall_File.iscell(:, 1) == 0, :) = [];

    if ii > 1
        % Behavior: 'Tif_N > 'CageID_Date_TrackN_sessionN.mat''
        behavior_fspec = 'floating_%d.mat';
        behavior_fname = sprintf(behavior_fspec,ii);
        floating = importdata(behavior_fname);
    end

    % Removing any artifacts in the data
    artifacts = find(isnan(floating.X));
    if sum(artifacts) > 0 % if there are artifacts
        floating.Y(artifacts) = [];
        floating.speed(artifacts) = [];
        floating.R(artifacts) = [];
        floating.alpha(artifacts) = [];
        data.DFF(:, artifacts) = [];
        data.spikes(:, artifacts) = [];
    end

    %% Functions
    %   Here we'll have the different functions you can pass your data through.
    %   Run whichever ones suite your data.
    % ----------------------------------------------------------------------- %
    % Parameters
    degree_size = 5;
    act_thresh = 0.03;
    sd_thresh = 2; % 3

    params = struct;
    params.oneD_bin_size = degree_size;
    params.act_thresh = act_thresh;
    params.sd_thresh = sd_thresh;
    % ----------------------------------------------------------------------- %

    % Functions
    [ data ] = Spike_Max(data);                                    % establishes maximum spike rate
    [ data ] = Normalizer(data);                                   % normalizes your activity
    [ data, floating ] = Moving_v3(1,data, floating);              % removes timepoints when mouse is not moving
    [ active_cells{ii,:} ] = Active_Cells(data, act_thresh);       % removes non-active cells from data
    [ data ] = DFF_transients(data, sd_thresh);

    DFF_trans{ii,:} = data.DFF_transients;

    % extract lap-by-lap activity and smooth traces
    [ activity_binned_1D_sde_smoothed, activity_binned_1D_smoothed{ii,:}, angle_bin_ids{ii,:}, ...
        consistency{ii,:}, lap_times{ii,:},  lap_by_lap_activity{ii,:} ] = ...
        OneD_Track_Anaysis_v2(data, floating, degree_size, data.trialID);

    % remove binned activity artifacts
    activity_binned_1D_sde_smoothed(isnan(activity_binned_1D_sde_smoothed)) = 0;
    activity_binned_1D_smoothed{ii,:}(isnan(activity_binned_1D_smoothed{ii,:})) = 0;
    
    % finding place cells using cohen's d
    [ place_cells, FWHM{ii,:}, params ] = PC_reliability_checker_WTR_v2(lap_by_lap_activity{ii,:}, 0, 0, 0, 0);
    
    [ spatial_info{ii,:}, ~ ] = SpatialInfoComputer(activity_binned_1D_smoothed{ii,:}, angle_bin_ids{ii,:}, lap_times{ii,:});

    [ speed_cells, speed_binned{ii,:}, spikes_speed_binned{ii,:}, sde_spikes_speed_binned{ii,:}, ~, ...
        ~, ~] = Speed_Cells(data, floating, lap_times{ii,:});

    PC_array{ii} = place_cells; % cell array of pc indices
    SC_array{ii} = speed_cells; % cell array of sc indices
end
data.DFF_transients = DFF_trans; % change to be transients from all environments

%% Get vector of all cells with a place field in at least one environment
% just take unique active cells, without repitition
all_pcs = [PC_array{1},PC_array{2},PC_array{3}];
valid_PCs = unique(all_pcs);
all_scs = [SC_array{1},SC_array{2},SC_array{3}];
valid_SCs = unique(all_scs);
all_nonactive = [find(~active_cells{1}),find(~active_cells{2}),...
    find(~active_cells{3})];
nonactive = unique(all_nonactive);
for n = 1:length(nonactive)
    curr_idx = nonactive(n);
    if ~isempty(find(valid_PCs == curr_idx)) % if the cell has been classified as a place cell
        valid_PCs(valid_PCs == curr_idx) = []; % remove non-active cells
    end
end


%% Plot according to input arguments
if plot_heatMap
    [valid_PCs,valid_SCs,sorted_PCs,PC_activity] = plotHeatMaps(PC_array,...
        SC_array,activity_binned_1D_smoothed,num_envs,valid_PCs,valid_SCs);
end
if plot_laps
   plotAllLapByLap(lap_by_lap_activity,valid_PCs);
end
if plot_oneD
    plotOneD(lap_by_lap_activity,valid_PCs,num_envs)
end
if plot_donut
    plotDonutHeatmap(valid_PCs,activity_binned_1D_smoothed);
end

%% Save data to current directory including date
if save_flag
    Save_Data(data,active_cells,spatial_info,activity_binned_1D_smoothed,...
        activity_binned_1D_sde_smoothed,lap_by_lap_activity,params,valid_PCs,...
        valid_SCs,FWHM,angle_bin_ids,lap_times,'_')
end
end