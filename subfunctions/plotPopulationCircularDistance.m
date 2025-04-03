function [scAngDist] = plotPopulationCircularDistance(combined,plot_flag)
%% Single Cell Angular Distance calculation: Figure 5C
%-------------------------------------------------------------------------%
%   Gets single cell angular distance across all recordings. Plots and runs
%   statistical model.
%
%   Written by NSW 04/05/2024 // Last updated by NSW 04/24/2024
%-------------------------------------------------------------------------%
if nargin < 1 || isempty(combined)
    disp('Select place cell responses...')
    [pname, ppath] = uigetfile('*.mat');
    combined = importdata(fullfile(ppath,pname));
end
if nargin < 2 || isempty(plot_flag)
    plot_flag = 1;
end

%% Set parameters
stages = {'D','P','E','M'};
n_envs = 3;
n_comps = 2;

%% Calculate angular distance
fun = @(x) x(1);
cellsz = cellfun(@size,combined,'uni',false);
all_cellnums = cellfun(fun,cellsz);
cellnums = sum(all_cellnums,2); % get number of every cell in array per stage
maxsz  = max(cellnums); % max number of cells in any stage

% initialize vectors for stats
stage_array = NaN(maxsz, size(combined,1)); % cells x stages

% initialize data for plotting
scAngDist = NaN(maxsz, n_comps, size(combined,1)); % cells x comps x stages
for ii = 1:size(combined,1) % for each stage
    stage_array(1:cellnums(ii),ii) = ones(cellnums(ii),1)*ii;
    curr_combined = cell2mat(combined(ii,:)'); % all cells for current stage
    AD_array = zeros(cellnums(ii),n_comps); % n cells in stage x n envs
    for dd = 1:cellnums(ii) % for each cell
        curr_rec = curr_combined(dd,:,:);
        curr_max = zeros(1,n_envs); % max peak for each environment
        for e = 1:n_envs
            curr_env = curr_rec(:,:,e);
            [~,max_vec] = max(curr_env, [], 2); % get vector of indices of max points
            curr_max(:,e) = max_vec;
        end
        r_score = abs(curr_max(:,1) - curr_max(:,2));
        s_score = abs(curr_max(:,1) - curr_max(:,3));

        % convert to circular distance
        if r_score > size(curr_env,2)/2
            r_score = abs(r_score - size(curr_env,2));
        end
        if s_score > size(curr_env,2)/2
            s_score = abs(s_score - size(curr_env,2));
        end

        AD_array(dd,1) = s_score;
        AD_array(dd,2) = r_score;
    end
    vector_comp = AD_array*5; % rescale to 1-180 degrees
    scAngDist(1:length(vector_comp),:,ii) = vector_comp;
end

%% Plot
if plot_flag
    figure
    for ii = 1:size(combined,1) % for each stage
        subplot(1,size(combined,1),ii);
        curr_stage = scAngDist(:,:,ii);
        groups = {'A>A''','A>B'};
        [p,~,~] = anova1(curr_stage,groups,'off');
        violinplot(curr_stage,groups);
        hold on
        ylim([0 180])
        ylabel('Angular Distance')
        title(stages{1,ii})
        subtitle(strcat('p = ', num2str(p)))
    end
    sgtitle('Single Cell Angular Distance')
end