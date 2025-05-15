function [ data, coords ] = Moving_mini2p(move_flag, move_thresh,...
    likelihood_thresh, data, coords, likelihood)
%-------------------------------------------------------------------------%
% MOVING_MINI2P Removes times where the mouse is not moving and coordinates 
% with low likelihood from the data.
%
%   [data, floating] = MOVING_V3(move_flag, data, floating) removes
%   times where the mouse is not moving from both the behavior data and
%   the suite2p activity data. This function works with both the newer
%   and older versions of the Neurotar software.
%
%   INPUT
%   'move_flag' is a boolean flag indicating whether to remove times
%   where the mouse is not moving. If 'move_flag' is not provided or
%   empty, it defaults to 1 (true).
%   'data' is a structure containing the suite2p activity data, with
%   fields 'DFF', 'DFF_raw', 'raw_F', and 'neuropil_F'.
%   'floating' is a structure containing the behavior data, with a field
%   'moving_times' that is a logical vector indicating the times when
%   the mouse is moving.
%
%   OUTPUT
%   The function returns 'data' and 'floating' with the times where the
%   mouse is not moving removed.
%
%   Written by NSW 11/16/2022 // Based on code from WTR Sept. 2020
%-------------------------------------------------------------------------%
if nargin < 1 || isempty(move_flag)
    move_flag = 1;
end
if nargin < 2 || isempty(move_thresh)
    move_thresh = 5;
end
if nargin < 3 || isempty(likelihood_thresh)
    likelihood_thresh = 0.3;
end

if move_flag
    %% Find indices where the mouse is not moving or likelihood is low
    % Compute Euclidean distance between consecutive coordinates
    diff_coords = diff(coords);
    movement = sqrt(sum(diff_coords.^2, 2));

    % Identify frames where movement is below threshold
    still_frames = movement < move_thresh;

    % Align with original coords (diff shortens by 1)
    still_frames = [false; still_frames];

    % Identify frames with high enough likelihood
    low_likelihood = likelihood <= likelihood_thresh;

    % Combine conditions: still and valid
    not_moving_time = still_frames & low_likelihood;

    %% Remove not moving time from activity data and coordinates
    coords(not_moving_time,:) = [];

    data.DFF(:, not_moving_time) = [];
    data.DFF_raw(:, not_moving_time) = [];
    data.spikes(:, not_moving_time) = [];
    data.raw_F(:, not_moving_time) = [];
    data.neuropil_F(:, not_moving_time) = [];
end