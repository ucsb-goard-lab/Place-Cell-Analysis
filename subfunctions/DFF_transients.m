function [data] = DFF_transients(data, sd_thresh)
%-------------------------------------------------------------------------%
% DFF_TRANSIENTS Finds the transients in the DFF of each neuron.
%
%   data = DFF_TRANSIENTS(data, sd_thresh) finds the transients in the
%   DFF (delta F/F) of each neuron in 'data' and masks all other points
%   to zero. This function follows the method used in Hainmuller and
%   Bartos 2018.
%
%   'data' is a structure with a field 'DFF', which is an NxM array
%   where N is the number of neurons and M is the number of time points.
%   'sd_thresh' is a scalar that specifies the threshold for defining
%   transients, in terms of standard deviations above the baseline.
%
%   The function returns 'data' with an additional field 'DFF_transients',
%   which is an NxM array of the same size as 'DFF'. This array contains
%   the DFF values at the times of the transients, and zeros elsewhere.
%
%   Written by WTR 12/27/20 // Last updated by WTR 01/10/21
%-------------------------------------------------------------------------%
N = size(data.DFF, 1);
DFF_transients = zeros(size(data.DFF));

for ii = 1:N
    sd = std(data.DFF(ii, :));
    baseline_times = find(data.DFF(ii, :) < (sd_thresh * sd));
    mean_baseline = mean(data.DFF(ii, baseline_times));
    sd_baseline = std(data.DFF(ii, baseline_times));
    transient_times = find(data.DFF(ii, :) > (sd_thresh * sd_baseline + mean_baseline));
    DFF_transients(ii, transient_times) = data.DFF(ii, transient_times);
end

data.DFF_transients = DFF_transients;