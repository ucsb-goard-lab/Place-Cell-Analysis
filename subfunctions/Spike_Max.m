function [ data ] = Spike_Max(data)
%-------------------------------------------------------------------------%
%SPIKE_MAX Implements a maximum spike rate cut-off as suggested to me by
%   Kevin Sit. 
%
%   Input:
%   data: a structure containing the spike data.
%
%   This function implements a maximum spike rate cut-off as suggested by
%   Kevin Sit. It caps the spikes of each cell at a maximum value determined
%   by the mean and standard deviation of the spikes. The maximum value is
%   set to be the mean plus three times the standard deviation. Any spikes
%   exceeding this value are set to the maximum value.
%
%   Output:
%   data: the input structure with the spikes capped at the maximum value.
%
%   Example:
%       data = Spike_Max(data)
%
%   See also SIZE, STD, MEAN.
%
%   Written by WTR 12/12/2019 // Last updated by WTR 12/19/2020
%-------------------------------------------------------------------------%
%% Globals
max_std = 3; 
N = size(data.spikes, 1); 

%% Capping the spikes
for ii = 1:N
    spikes_ii = data.spikes(ii, :); 
    std_ii = std(spikes_ii);
    mean_ii = mean(spikes_ii);
    spikes_ii(spikes_ii > (mean_ii + std_ii * max_std)) = mean_ii + std_ii * max_std;
    data.spikes(ii, :) = spikes_ii; 
end

