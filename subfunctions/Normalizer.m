function [ data ] = Normalizer(data)
%-------------------------------------------------------------------------%
%   This script normalizes the data. It will either normalize the data
%   against itself, or another recording (so that single recording sessions
%   can be appropriately compared). 
%
%    'data' is a structure with fields 'spikes' and 'DFF', which are arrays
%   containing the spike and DFF (delta F/F) data, respectively
%   Written by WTR 09/19/2019 // Last updated by WTR 3/15/2020
%-------------------------------------------------------------------------%
%% Normalizing 
data.spikes = data.spikes ./ max(data.spikes, [], 2); 
data.DFF = data.DFF ./ max(data.DFF, [], 2);

end

