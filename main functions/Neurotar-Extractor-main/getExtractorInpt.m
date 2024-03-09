function [frames, short_frames, rrate, fname_bank]  = getExtractorInpt()
%GETEXTRACTORINPT Extracts necessary input parameters for a NeurotarExtractor.
%   This function calculates the number of frames from 2p recordings and
%   the number of frames that correspond to the Neurotar recording. It also
%   extracts the exact frame rate from the .env file.
%
%   Outputs:
%   - frames: A vector containing the number of frames for each image.
%   - short_frames: A vector containing the number of frames that correspond to the Neurotar recording.
%   - rrate: A vector containing the exact frame rate for each recording.
%   - fname_bank: A cell array containing the filenames of the .mat files.

fname_bank = {'floating_1.mat', 'floating_2.mat', 'floating_3.mat'}; % specify filenames

% get raw number of frames from 2p recording
tifs = dir('*.tif');
if ~isempty(tifs)
    disp('Select your tif files')
    tif_files = uigetfile('*.tif','MultiSelect','on');
    num_files = length(tif_files);
    disp('Calculating recording frames...')
    frames = zeros(1,num_files);
    names = cell(1,num_files);
    for ii = 1:num_files
        info = imfinfo(tif_files{ii});
        numberOfPages = length(info);
        names{1,ii} = info(ii).Filename;
        frames(1,ii) = numberOfPages; % calculate number of frames for each image
    end
else % for later processing steps where tifs have been deleted
    disp('No tifs found: extracting frames from Fall files...')
    fall_files = dir('Fall_*');
    num_files = length(fall_files);
    frames = zeros(1,num_files);
    for ii = 1:num_files
        Fall = importdata(fall_files(ii).name);
        frames(1,ii) = size(Fall.F,2);
    end
end
tifdir = pwd;
cd ..

% get the number of frames that correspond to the neurotar recording
short_frames = zeros(1,num_files);
rrate = zeros(1,num_files);
xml_dirs = dir('TSeries*');
for ee = 1:num_files
    curr_xmldir = xml_dirs(ee).name;
    cd(curr_xmldir)
    f = NewNeurotarExtractor([], isMoving = true, RECORDING_FRAMES = frames(ee),...
        RECORDING_RATE = 10);
    % f = importdata(strcat('AllTifs\',fname_bank{ee}));
    ndata = f.data;
    HW_timestamp = ndata.HW_timestamp; % extract relative timestamps
    last_time = double(HW_timestamp(end))/1000; % get end time
    if last_time == 0 % sometimes last frame gets rewritten to 0
        last_time = double(HW_timestamp(end-1))/1000; % take time right before end
    end

    % get exact frame number from xml file
    xmlname = dir('*.xml'); % select xml file
    xml = readstruct(xmlname.name);
    relative_time = zeros(length(xml.Sequence.Frame),1);
    for rr = 1:length(xml.Sequence.Frame)
        relative_time(rr,1) = xml.Sequence.Frame(rr).relativeTimeAttribute;
    end
    [~,cutoff]=min(abs(relative_time-last_time)); % get closest frame to cutoff point
    if cutoff < 500 % if there is an error with the xml file
        cutoff = round(last_time*10);
        disp('Error in xml file: extracing cutoff time from tdms')
    end
    short_frames(ee) = cutoff;

    % get exact frame rate from env file
    envname = dir('*.env'); % select xml file
    env_file = importdata(envname.name);
    for j = 1:length(env_file)
        if strfind(env_file{j},'repetitionPeriod') > 0
            cfg_line = env_file{j};
            index = strfind(cfg_line,'repetitionPeriod');
            frameRate = 1/sscanf(cfg_line(index:end),'repetitionPeriod="%f"');
        end
    end
    rrate(ee) = frameRate;

    cd ..
end
cd(tifdir)