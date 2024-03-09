function [data] = spikeInference(filename,save_flag)
%SPIKEINFERENCE Infers spike rates from calcium imaging data.
%
%   inputs:
%   filename: a string containing the name of the .mat file to load.
%   save_flag: a binary flag indicating whether to save the data (1) or not (0).
%
%   This function infers spike rates from calcium imaging data using the
%   deconvolveCa function. If no inputs are provided, the function prompts
%   the user to select a .mat file to load and sets the save_flag to 1. If
%   only the filename is provided, the function sets the save_flag to 1.
%
%   The function loads the data from the .mat file, infers the spikes for
%   each cell in the data, and adds the inferred spikes to the data
%   structure. If save_flag is set to 1, the function saves the data
%   structure to the .mat file.
%
%   output:
%   data: the input structure with the inferred spikes added.
%
%   Example:
%       data = spikeInference('data.mat', 1)
%
%   See also UIGETFILE, LOAD, DISP, ZEROS, SIZE, DECONVOLVECA, SAVE.
%
% Updated 09Jul2019 changed the save from eval(['save' filename 'data']) to current, fixing bad form

if nargin == 0
    [filename,pathname] = uigetfile('.mat');
    load(filename);
    save_flag = 1;
elseif nargin == 1
    save_flag = 1;
end

load(filename);
disp('Inferring spikes...')
dffDeconv = zeros(size(data.DFF));

for n = 1:size(data.DFF, 1)
    % get trace and run deconvolution
    trace = data.DFF(n,:);
    [denoised,spikes,opt] = deconvolveCa(trace, 'ar1' ,'foopsi', 'optimize_pars');
    dffDeconv(n,:) = spikes;
end
data.spikes = dffDeconv;
close all

if save_flag
save(filename,'data');
end
end

