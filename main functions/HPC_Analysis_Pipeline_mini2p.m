function [data, activityMap, valid_PCs] = ...
    HPC_Analysis_Pipeline_mini2p(data, Fall, dlc, save_flag, plot_flag)
% HPC Analysis Pipeline for mini2p 
%-------------------------------------------------------------------------%
%   This function is an updated version of the HPC analysis pipeline that
%   uses Cohen's d for place cell (PC) thresholding. Should be run from the
%   directory containing GOARD_method_processed_data.mat, the mp4 video, 
%   and the 'suite2p' analyzed data folder.
%
%   Inputs
%   - data: The processed data output from suite2p.
%   - dlc: Deeplabcut tracking file, here use the csv output.
%   - save_flag: A flag to save data.
%   - plot_flag: A flag to plot data.
%
%   Outputs:
%   - data: Processed suite2p data.
%   - activityMap: A 3D array containing binned cell DFF.
%   - valid_PCs: A vector containing the indices of all place cells.
%
%   Written by NSW 05/24/2022 // Last updated by NSW 05/15/2025
%-------------------------------------------------------------------------%
if nargin < 1 || isempty(data)
    data = importdata('GOARD_method_processed_data.mat');
end
if nargin < 2 || isempty(dlc)
    dname = dir('*.csv');
    dlc = importdata(dname(1).name); 
end
if nargin < 3 || isempty(Fall)
    Fall = importdata('suite2p/plane0/Fall.mat');
end
if nargin < 4 || isempty(save_flag)
    save_flag = 1; 
end
if nargin < 5 || isempty(plot_flag)
    plot_flag = 1;
end

%% Delete DFF traces and spikes from non-cells
data.DFF(Fall.iscell(:, 1) == 0, :) = []; 
data.spikes(Fall.iscell(:, 1) == 0, :) = [];
data.DFF_raw(Fall.iscell(:, 1) == 0, :) = [];

%% Extract relevant variables from DLC data
coords = dlc.data(:,2:3); % nose coordinates
likelihood = dlc.data(:,4); % likelihood

%% Functions
%   Here we'll have the different functions you can pass your data through.
%   Run whichever ones suite your data.
% ----------------------------------------------------------------------- %
% Parameters
bin_size = 131; % N x N bin size in pixels. 131 = 3cm
sd_thresh = 2; % 3
move_thresh = 5; % number of pixels between coordinates below which is considered 'not moving'
likelihood_thresh = 0.3; % dlc likelihood below which to remove coordinates

% ----------------------------------------------------------------------- %

% Functions
[ data ] = Spike_Max(data);                              % establishes maximum spike rate
[ data ] = Normalizer(data);                             % normalizes your activity
[ data, coords ] = Moving_mini2p(1, move_thresh,...
    likelihood_thresh, data, coords, likelihood);        % removes timepoints when mouse is not moving or when dlc tracking isn't good
[ data ] = DFF_transients(data, sd_thresh);

% extract lap-by-lap activity and smooth traces
[ activityMap ] = binActivity(data, coords, bin_size);

% finding place cells using cohen's d
[ valid_PCs ] = PC_reliability_checker_mini2p(data, bin_size, coords);

%% Plot according to input arguments
if plot_flag
    plotPCHeatmaps(activityMap, valid_PCs)
end

%% Save data to current directory including date
if save_flag
    Save_Data_mini2p(data, activityMap, valid_PCs)
end
