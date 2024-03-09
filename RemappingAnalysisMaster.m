%%% Master script for running place cell remapping pipeline

%% Add all functions to path
basedir = fileparts(which('RemappingAnalysisMaster.m')); 
addpath(genpath(basedir));

%% Extract data 
envs = 3; % number of environments
[frames, short_frames, rrate, fname_bank] = getExtractorInpt(); 
for ii = 1:envs
    Neurotar = NewNeurotarExtractor([], isMoving = true, RECORDING_FRAMES = short_frames(ii),...
        RECORDING_RATE = rrate(ii)); % extract neurotar data from tdms file, inpt 1 = rec name
    save(fname_bank{ii},'Neurotar') % save under filename defined in getExtractorInpt
end
data = suite2p2data([],[],1,1,envs,frames,short_frames); % extract from suite2p; inpts = (F, Fneu, save_flag, deconcat_flag, num_envs, frames)

%% Analyze place cells
% inputs: (num_envs, data, floating, plot_donut, plot_laps, plot_oneD, plot_heatMap, save_flag)
HPC_Analysis_Pipeline_Method3(envs, [], [], 0, 0, 0, 1, 1);

