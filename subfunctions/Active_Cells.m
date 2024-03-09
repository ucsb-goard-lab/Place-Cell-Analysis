function [ active_cells ] = Active_Cells(data, act_thresh)
%-------------------------------------------------------------------------%
%ACTIVE_CELLS Computes the number of active cells in the data.
%   This function computes the number of active cells in the data, following the method used in 
%   Hainmueller and Bartos 2018. It calculates a threshold for each cell based on its mean and standard 
%   deviation of DF/F values. It then finds the transients that exceed this threshold. If the number of 
%   transients exceeds a certain threshold (act_thresh), the cell is considered active.
%
%   Inputs:
%   - data: A structure containing the processed data. It includes the DF/F traces.
%   - act_thresh: The threshold for the number of transients to consider a cell as active.
%
%   Output:
%   - active_cells: A vector containing 1 for active cells and 0 for non-active cells.
%
%   Written by WTR 08/14/2020 // Last updated by NSW 07/31/2023
%-------------------------------------------------------------------------%
%%
active_cells = zeros(1, size(data.DFF, 1)); 
sigma = 2;

%%
for ii = 1:size(data.DFF, 1)
    transient_thresh = mean(data.DFF(ii, :)) + sigma * std(data.DFF(ii, :)); 
    transients = find(data.DFF(ii, :) > transient_thresh); 
    
    if length(transients)  / 10 > act_thresh
        active_cells(ii) = 1;
    end
end