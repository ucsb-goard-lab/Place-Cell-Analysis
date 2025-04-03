function [pvecCorr_bootstrap_vecs] = plotPVec_Population_Bootstrap(iterations,combinedData,plot_flag)
%% PVec calculation: Figure 5D
%-------------------------------------------------------------------------%
%   Uploads one-dimensional activity by recording and gets population
%   vectors, acconting for low-cell count recordings via boostrapping.
%
%   Written by NSW 04/24/2024 // Last updated by NSW 04/24/2024
%-------------------------------------------------------------------------%
if nargin < 1 || isempty(iterations)
    iterations = 100;
end
if nargin < 2 || isempty(combinedData)
    disp('Select place cell responses...')
    [pname, ppath] = uigetfile('*.mat');
    combinedData = importdata(fullfile(ppath,pname));
end
if nargin < 3 || isempty(plot_flag)
    plot_flag = 1;
end

% Run bootstrapping: End result is n correlation averages for each stage in
% each iteration, where n is the # of recordings in that stage.
% Correlations are taken across columns of responses
% First row is stability corr (A > A'), second row is remapping corr (A > B)
n_comps = 2; % number of comparisons (A > A'and A > B)

%% extract matrix with number of cells per mouse/stage, calculate minimum number of cells to use
numCelMat = zeros(size(combinedData,1),size(combinedData,2)); % stages x recs
for s = 1:size(combinedData,1) % for each stage
    for m = 1:size(combinedData,2) % for each recording
        numCelMat(s,m) = size(combinedData{s,m},1);
    end
end
numCell = min(numCelMat);        % mininum number of cells across stages per mouse (to keep stages equal)

pvecCorr_bootstrap = zeros(size(combinedData,2),n_comps,size(combinedData,1),iterations); % recs x comparisons x stages x iterations
pvecCorr_bootstrap_vecs = zeros(size(combinedData,2)*iterations,n_comps,size(combinedData,1)); % vectorized version for plotting

stab_CI_bootstrap = zeros(size(combinedData,2)*iterations,n_comps,size(combinedData,1));
rem_CI_bootstrap = zeros(size(combinedData,2)*iterations,n_comps,size(combinedData,1));
count = 1;
for i = 1:iterations
    pvecCorr_all = zeros(size(combinedData,2),n_comps,size(combinedData,1)); % recs x comparisons x stages
    stab_CI_all = zeros(size(combinedData,2),n_comps,size(combinedData,1));
    rem_CI_all = zeros(size(combinedData,2),n_comps,size(combinedData,1));
    % pvecCorr = zeros(n_comps,size(combinedData{1},2),length(combinedData)); % dims: n_comps x n_bins x n_stages
    for ii = 1:size(combinedData,1) % for each stage
        curr_stage = combinedData(ii,:); % all recordings in current stage
        corrs = zeros(size(curr_stage,2),n_comps); % final correlations are number of recordings x n_comps
        stab_CI = zeros(size(curr_stage,2),n_comps); % confidence intervals for stability: lower bound, upper bound
        rem_CI = zeros(size(curr_stage,2),n_comps); % confidence intervals for remapping: lower bound, upper bound
        for aa = 1:size(curr_stage,2) % for each recording within a stage
            curr_array = curr_stage{aa};
            curr_numcell = numCell(aa);
            rand_cells = randperm(size(curr_array,1),curr_numcell); % X unique integers selected from cell ID's
            curr_array = curr_array(rand_cells,:,:); % just take randomly selected subset
            pvecCorr = zeros(n_comps,size(curr_array,2)); % 2 comparisons x 72 bins
            for cc = 1:n_comps
                for j = 1:size(curr_array,2) % 1 to num bins
                    vec1 = curr_array(:,j,1); % population vector at current bin for env A
                    vec2 = curr_array(:,j,2+(abs(cc-n_comps))); % population vector at current bin for env A' for i == 1 and B for i == 2
                    R = corrcoef(vec1,vec2); % get pearson's correlation
                    pvecCorr(cc,j) = R(1,2);
                end
            end
            stab_ste = std(pvecCorr(1,:))/sqrt(length(pvecCorr(1,:))); % get standard error
            stab_ts = tinv([0.025  0.975],length(pvecCorr(1,:))-1); % get T-Score
            stab_CI(aa,:) = mean(pvecCorr(1,:)) + stab_ts*stab_ste; % calculate Confidence Intervals

            rem_ste = std(pvecCorr(2,:))/sqrt(length(pvecCorr(2,:))); % get standard error
            rem_ts = tinv([0.025  0.975],length(pvecCorr(2,:))-1); % get T-Score
            rem_CI(aa,:) = mean(pvecCorr(2,:)) + rem_ts*rem_ste; % calculate Confidence Intervals

            corrs(aa,:) = mean(pvecCorr,2,'omitnan')';
        end

        pvecCorr_all(:,:,ii) = corrs;
        stab_CI_all(:,:,ii) = stab_CI;
        rem_CI_all(:,:,ii) = rem_CI;
    end
    pvecCorr_bootstrap(:,:,:,i) = pvecCorr_all;
    pvecCorr_bootstrap_vecs(count:count+(size(curr_stage,2)-1),:,:) = pvecCorr_all;

    stab_CI_bootstrap(count:count+(size(curr_stage,2)-1),:,:) = stab_CI_all;
    rem_CI_bootstrap(count:count+(size(curr_stage,2)-1),:,:) = rem_CI_all;
    count = count + size(curr_stage,2);

    if rem(i,10)==0 % display current iteration
        disp(['iteration ' num2str(i) '/' num2str(iterations) ' complete'])
    end   
end

if plot_flag
    % plot subplots comparring PVec corrs in violin plots
    addpath('E:\Code\Plotting')
    stages = {'D','P','E','M'};
    groups = {'A>A''','A>B'};
    colors = {[0.05 0.07 0.57],[1,0.38,0.18]};
    
    figure
    for ii = 1:size(pvecCorr_bootstrap,2)%size(combinedData,1) % for each comparison
        subplot(1,size(pvecCorr_bootstrap,2),ii);
        % curr_plot = pvecCorr_bootstrap_vecs(:,:,ii); % first column is stability, second is remapping
        % curr_plot = pvecCorr_byRecording(:,:,ii); 
        curr_plot = squeeze(pvecCorr_bootstrap_vecs(:,ii,:));
        if ii == 1 % for stability
            all_CI = stab_CI_bootstrap;
        else % for remapping
            all_CI = rem_CI_bootstrap;
        end
        mean_CI = squeeze(mean(all_CI,'omitmissing')); % lower and upper bound 95% CI for all stages
        
        for ss = 1:size(pvecCorr_bootstrap,3) % for all stages
            [density, value] = ksdensity(curr_plot(:,ss));
            density = density(value >= min(curr_plot(:,ss)) & value <= max(curr_plot(:,ss)));
            value = value(value >= min(curr_plot(:,ss)) & value <= max(curr_plot(:,ss)));
            value(1) = min(curr_plot(:,ss));
            value(end) = max(curr_plot(:,ss));

            width = 0.3/max(density);

            % plot the violin
            fill([ss+density*width ss-density(end:-1:1)*width], ...
                [value value(end:-1:1)], colors{ii},...
                'EdgeColor', colors{ii});

            % plot upper and lower CI
            quartiles = [mean_CI(1,ss),median(curr_plot(:,ss)),mean_CI(2,ss)];
            hold on
            fill(ss+[-1,1,1,-1]*.02, ...
                [quartiles(1) quartiles(1) quartiles(3) quartiles(3)],...
                [0 0 0], 'EdgeColor', [0 0 0]);
            scatter(ss,quartiles(2),30,[1 1 1],'filled','black') % plot mean
            scatter(ss,quartiles(2),20,[1 1 1],'filled','white') % plot mean
        end

        hold on
        ylim([-0.5 0.6])
        xlim([0.5 4.5])
        ylabel('Spatial Correlation')
        title(groups{ii})
        yline(0,'--')
    end
    sgtitle('Spatial Population Vector Correlation')
end

end