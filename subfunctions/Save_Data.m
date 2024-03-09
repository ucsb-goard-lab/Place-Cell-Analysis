function Save_Data(data,active_cells,spatial_info,activity_binned_1D_smoothed,...
    activity_binned_1D_sde_smoothed,lap_by_lap_activity,params,valid_PCs,...
    valid_SCs,FWHM,angle_bin_ids,lap_times,fileext)
%SAVE_DATA Saves processed data to a .mat file.
%
%   Takes twelve inputs:
%   data: a structure containing the processed data.
%   active_cells: a vector containing the indices of active cells.
%   spatial_info: a structure containing spatial information.
%   activity_binned_1D_smoothed: a 3D array containing smoothed 1D activity.
%   activity_binned_1D_sde_smoothed: a 3D array containing smoothed 1D activity standard error.
%   lap_by_lap_activity: a cell array where each cell contains a 3D array of lap by lap activity.
%   params: a structure containing parameters.
%   valid_PCs: a vector containing the indices of valid place cells.
%   valid_SCs: a vector containing the indices of valid speed cells.
%   FWHM: a vector containing the full width at half maximum for each cell.
%   angle_bin_ids: a vector containing the bin IDs for each angle.
%   lap_times: a vector containing the lap times.
%   fileext: a string containing the file extension.
%
%   This function saves the processed data to a .mat file. The file is saved
%   to the current directory with a filename that includes the current date
%   (in the format YYMMDD), the method used to process the data (either
%   'suite2p' or 'GoardMethod'), and the given file extension.
%
%   If no data is provided, the function will prompt the user to select a
%   .mat file to import. If no file extension is provided, the function
%   will set the file extension to '_'.
%
%   Example:
%       Save_Data(data,active_cells,spatial_info,activity_binned_1D_smoothed,...
%       activity_binned_1D_sde_smoothed,lap_by_lap_activity,params,valid_PCs,...
%       valid_SCs,FWHM,angle_bin_ids,lap_times,'_')
%
%   See also UIGETFILE, IMPORTDATA, ISFIELD, DATETIME, YMD, SAVE.

 if nargin < 1 || isempty(data)
     disp('Select processed data')
     d_fname = uigetfile('*.mat');
     data = importdata(d_fname);
 end
 if nargin < 12 || isempty(fileext)
     fileext = '_'; % set file extension to blank if none is given
 end

 % if length(active_cells) == 3 % for remapping experiments
 %     data.isCells = sum(active_cells{1});
 % else
 %     data.isCells = sum(active_cells);
 % end
data.spatial_info = spatial_info;
data.oneD_activity = activity_binned_1D_smoothed;
data.oneD_activity_SEM = activity_binned_1D_sde_smoothed;
data.lap_by_lap_activity = lap_by_lap_activity;
data.params = params;
data.valid_pcs = valid_PCs;
data.active_cells = active_cells;
% data.sorted_pcs = sorted_PCs;
% data.PC_activity = PC_activity;
data.valid_SCs = valid_SCs;
data.FWHM = FWHM;
data.binned_angles = angle_bin_ids;
data.lap_times = lap_times;
% data.speed_score = speed_score;
% data.speed_binned = speed_binned;
% data.speed_activity = spikes_speed_binned;
% data.speed_activity_SEM = sde_spikes_speed_binned;
% data.n_thresh = n_thresh;
% data.p_thresh = p_thresh;
% data.angles = angle_bin_ids;
% data.consistency = consistency;
% data.mean_rate = mean_rate;
% data.PCs = PCs;
% data.speed_cells = SCs;

[y,m,d] = ymd(datetime("today"));
y = num2str(y);
m = num2str(m);
d = num2str(d); 
if length(d) == 1 % if it's a one-digit day
    d = strcat('0',d);
end
if length(m) == 1
    m = strcat('0',m);
end
if isfield(data,'spikes') % because only suite2p data outputs spikes
    methodflag = 'suite2p';
else
    methodflag = 'GoardMethod';
end
save(strcat(y(3:end),m,d,'_processed_data_Cohens_1.2_',methodflag,fileext,'NW','.mat'), 'data'); % saves to file with current date (YYMMDD)
disp('Done! Processed data saved to current directory.')
end