function [ spatial_info, mean_rate ] = SpatialInfoComputer(activity, angles, lap_times)
%-------------------------------------------------------------------------%
%SPATIALINFOCOMPUTER Computes the spatial information and mean rate.
%
%   Inputs:
%   activity: a 2D array containing the activity data.
%   angles: a vector containing the angles.
%   lap_times: a 2D array containing the lap times.
%
%   This function first normalizes the activity data. It then computes the
%   counts of each angle within the lap times. Finally, it computes the
%   spatial information and mean rate for each row of the activity data
%   using the Spatial_Information_v2 function.
%
%   Outputs:
%   spatial_info: a vector containing the computed spatial information for each row of the activity data.
%   mean_rate: a vector containing the computed mean rate for each row of the activity data.
%
%   Example:
%       [spatial_info, mean_rate] = SpatialInfoComputer(activity, angles, lap_times)
%
%   See also MAX, SUM, SPATIAL_INFORMATION_V2.
%
%   Written by WTR 02/06/2021 // Last updated by WTR 02/06/2021
%-------------------------------------------------------------------------%
%% Normalizing activity 
activity = activity ./ max(activity, [], 2); 

%%
angles2 = [];
for nn = 1:size(lap_times, 1)
    angles2 = [angles2, angles(lap_times(nn, 1):lap_times(nn, 2))];
end 

counts = zeros(1, size(activity, 2)); 
for ii = 1:size(activity, 2)
    counts(ii) = sum(angles2 == ii); 
end

spatial_info = zeros(1, size(activity, 1)); 
mean_rate = zeros(1, size(activity, 1)); 

for ii = 1:size(activity, 1)
    [spatial_info(ii), mean_rate(ii)] = Spatial_Information_v2(activity(ii, :), counts);
end



