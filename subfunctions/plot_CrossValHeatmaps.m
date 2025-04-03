function [] = plot_CrossValHeatmaps(lap_by_lap)
%% Plot cross-validated heatmaps: Figure 5B
%-------------------------------------------------------------------------%
%   Gets lap-by-lap responses for place cells from env. A , cross-validates
%   them, sorts by max peak in env. A, and plots as heat maps for each
%   environment in each stage.
%
%   Written by MJG 02/01/2024 // Last updated by NSW 04/25/2024
%-------------------------------------------------------------------------%
if nargin < 1 || isempty(lap_by_lap)
    disp('Select place cell responses...')
    [lname, lpath] = uigetfile('*.mat');
    lap_by_lap = importdata(fullfile(lpath,lname));
end

% Set parameters
num_stages = 4;
num_envs = 3;

%% Initialize matrices for concatenated data
resp_mat_A_allcells = cell(1,num_stages);
resp_mat_B_allcells = cell(1,num_stages);
resp_mat_Ap_allcells = cell(1,num_stages);
pref_idx_allcells = cell(1,num_stages);
pref_idx_byRecording = cell(size(lap_by_lap));
for aa = 1:size(lap_by_lap,2) % for each animal
    curr_animal = lap_by_lap(:,aa);
    for ss = 1:length(curr_animal) % for each stage
        curr_stage = curr_animal{ss};

        %% Extract parameters
        lap = size(curr_stage{1},1); % # laps in env A
        pos = size(curr_stage{1},2); % # bins in env A
        cel = size(curr_stage{1},3); % # cells in env A

        %% Determine preferred position using odd trials and calculate sorting vec
        for j = 1:cel
            resp_A = curr_stage{1}(:,:,j);
            mean_odd_resp = mean(resp_A([1:2:lap],:),1);
            [~,pref_idx(j)] = max(mean_odd_resp);
        end

        %% Cross-validate env A response using even trials, avg trials for env B/Ap
        resp_mat_A = squeeze(mean(curr_stage{1}([2:2:lap],:,:),1))';
        resp_mat_B = squeeze(mean(curr_stage{2},1))';
        resp_mat_Ap = squeeze(mean(curr_stage{3},1))';

        %% Add data to concatenated matrices
        resp_mat_A_allcells{1,ss} = [resp_mat_A_allcells{1,ss}; resp_mat_A];
        resp_mat_B_allcells{1,ss} = [resp_mat_B_allcells{1,ss}; resp_mat_B];
        resp_mat_Ap_allcells{1,ss} = [resp_mat_Ap_allcells{1,ss}; resp_mat_Ap];
        pref_idx_allcells{1,ss} = [pref_idx_allcells{1,ss}, pref_idx];
        pref_idx_byRecording{ss,aa} = pref_idx;

        clearvars resp_mat_A resp_mat_B resp_mat_Ap pref_idx
    end
end

%% determine sorting index for concatenated cells
sortIdx = cell(size(pref_idx_allcells));
for si = 1:num_stages
    curr_pref = pref_idx_allcells{si};
    [~,sortIdx{si}] = sort(curr_pref);
end

%% Sort by latency in env A preferred position
sort_mat_A = cell(size(sortIdx));
sort_mat_B = cell(size(sortIdx));
sort_mat_Ap = cell(size(sortIdx));
for sm = 1:num_stages
    curr_sort = sortIdx{sm};

    curr_respA = resp_mat_A_allcells{sm};
    sort_mat_A{sm} = curr_respA(curr_sort,:);

    curr_respB = resp_mat_B_allcells{sm};
    sort_mat_B{sm} = curr_respB(curr_sort,:);

    curr_respAp = resp_mat_Ap_allcells{sm};
    sort_mat_Ap{sm} = curr_respAp(curr_sort,:);
end

%% Plot
figure
all_sorted = [sort_mat_A;sort_mat_B;sort_mat_Ap];
titles = [{'Environment A'},{'Environment B'},{'Environment Ap'}];
colormats = [{[linspace(1,0,100)' linspace(1,0,100)' linspace(1,0.6,100)']},
        {[linspace(1,0.7,100)' linspace(1,0.2,100)' linspace(1,0,100)']},
        {[linspace(1,0,100)' linspace(1,0,100)' linspace(1,0.3,100)']}];
count = 0;
for ii = 1:num_stages
    for i = 1:num_envs
        curr_sorted = all_sorted{i,ii}; % load
        curr_sorted(curr_sorted>0.4) = max(max(curr_sorted)); % increase contrast
        curr_sorted = rescale(curr_sorted); % normalize
        count = count+1;
        ax(count) = subplot(num_stages,3,count);
        imagesc(curr_sorted)
        title(titles{i})
        ylabel('Place cell # (sorted)')
        xlabel('Position (cm)')

        colormat = colormats{i}; 
        colormap(ax(count),colormat)
        colorbar
        set(gcf,'color',[1 1 1])
    end
end