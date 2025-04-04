function [non_overlapping_CI_ap,non_overlapping_CI_b] = PC_Position_Decoder_Bootstrap()
%%% Run multipe iterations of PC_Position_Decoder_allMice
%%% to calculate average errors across stages
%   Written: MG 231005, Latest update: NW 240611

%% find data
disp('Select place cell responses...')
[lname, lpath] = uigetfile('*.mat');
fname = fullfile(lpath,lname);

%% params
iterations = 100;    % number of times to iterate
plotFlag = 0;        % turn off plotting to reduce processing time

stageList = {'Diestrus','Proestrus','Estrus','Metestrus'};
numStage = length(stageList);

disp('Responses selected, beginning decoding.')
[plotResults] = PC_Position_Decoder_allMice(plotFlag, fname);
numPos = plotResults.numPos;


%% initialize
error_B_avg = zeros(numPos,numStage,iterations);
error_Ap_avg = zeros(numPos,numStage,iterations);
likelihoodMat_B = zeros(numPos,numPos,numStage,iterations);
likelihoodMat_Ap = zeros(numPos,numPos,numStage,iterations);
CI_mat_Ap = zeros(2,numStage,iterations); % 2 for upper and lower bound
CI_mat_B = zeros(2,numStage,iterations);


%% run iterations
for i = 1:iterations

    [plotResults] = PC_Position_Decoder_allMice(plotFlag, fname);

    error_B_avg(:,:,i) = plotResults.error_B;
    error_Ap_avg(:,:,i) = plotResults.error_Ap;
    likelihoodMat_B(:,:,:,i) = plotResults.likelihood_B;
    likelihoodMat_Ap(:,:,:,i) = plotResults.likelihood_Ap;

    for ii = 1:numStage
        ap_error = plotResults.error_Ap(:,ii);
        b_error = plotResults.error_B(:,ii);

        Ap_ste = std(ap_error)/sqrt(length(ap_error)); % get standard error
        Ap_ts = tinv([0.025  0.975],length(ap_error)-1); % get T-Score
        CI_mat_Ap(:,ii,i) = mean(ap_error) + Ap_ts*Ap_ste; % calculate Confidence Intervals

        B_ste = std(b_error)/sqrt(length(b_error)); % get standard error
        B_ts = tinv([0.025  0.975],length(b_error)-1); % get T-Score
        CI_mat_B(:,ii,i) = mean(b_error) + B_ts*B_ste; % calculate Confidence Intervals
    end

    if rem(i,10)==0
        disp(['iteration ' num2str(i) '/' num2str(iterations) ' complete'])
    end
end
CI_arrays = cat(4,CI_mat_Ap,CI_mat_B);


%% calculate mean likelihood and bootstrap distributions
likelihoodMat_B_mean = mean(likelihoodMat_B,4);
likelihoodMat_Ap_mean = mean(likelihoodMat_Ap,4);
error_bootstrapped = cat(3,squeeze(mean(error_Ap_avg,1)),...
    squeeze(mean(error_B_avg,1))); % 3rd dim: Ap > B

%% normalize probability
min_prob = min([min(min(min(likelihoodMat_B_mean))) min(min(min(likelihoodMat_Ap_mean)))]);
likelihoodMat_B_norm = likelihoodMat_B_mean-min_prob;
likelihoodMat_Ap_norm = likelihoodMat_Ap_mean-min_prob;
max_prob = max([max(max(max(likelihoodMat_B_norm))) max(max(max(likelihoodMat_Ap_norm)))]);
likelihoodMat_B_norm = likelihoodMat_B_norm/max_prob;
likelihoodMat_Ap_norm = likelihoodMat_Ap_norm/max_prob;
new_max_prob = max([max(max(max(likelihoodMat_B_norm))) max(max(max(likelihoodMat_Ap_norm)))]);
new_min_prob = min([min(min(min(likelihoodMat_B_norm))) min(min(min(likelihoodMat_Ap_norm)))]);

%% plot mean likelihood matrices

figure(1)

for s = 1:numStage

    subplot(2,4,s);
    imagesc(likelihoodMat_Ap_norm(:,:,s))
    % cmap = colormapMaker([255,255,255;58,74,159]);
    colormap parula
    colorbar
    % clim([s_min s_max])
    clim([new_min_prob new_max_prob])
    set(gca,'YDir','normal')
    hold on
    plot(1:numPos,'k:')
    xlabel('Actual position')
    ylabel('Estimated position')% (maximum likelihood)')
    set(gca,'FontSize',14);
    axis square
    title(['A>A_p_r: ' stageList{s}])

    subplot(2,4,s+4);
    imagesc(likelihoodMat_B_norm(:,:,s))
    % cmap = colormapMaker([255,255,255;245,127,67]);
    colormap parula
    colorbar
    % clim([r_min r_max])
    clim([new_min_prob new_max_prob])
    set(gca,'YDir','normal')
    hold on
    plot(1:numPos,'k:')
    xlabel('Actual position')
    ylabel('Estimated position')% (maximum likelihood)')
    title(['A>B: ' stageList{s}])% ' stage'])
    set(gca,'FontSize',14);
    axis square

end

set(gcf,'Position',[150 150 1250 550])
set(gcf,'color',[1 1 1])


%% Plot error distributions

figure(2)

categories = [{'A>A'''},{'A>B'}];
colors = {[0.05 0.07 0.57],[1,0.38,0.18]};
for s = 1:length(categories) % Ap then B

    subplot(1,2,s)
    curr_env_error = error_bootstrapped(:,:,s);
    curr_env_CI = CI_arrays(:,:,:,s);
    mean_CI = mean(curr_env_CI,3); % lower bound; upper bound, for each stage

    for ss = 1:numStage
        curr_error = curr_env_error(ss,:)';
        curr_CI = mean_CI(:,ss); % upper bound; lower bound

        [density, value] = ksdensity(curr_error);
        density = density(value >= min(curr_error) & value <= max(curr_error));
        value = value(value >= min(curr_error) & value <= max(curr_error));
        value(1) = min(curr_error);
        value(end) = max(curr_error);

        width = 0.3/max(density);

        % plot the violin
        fill([ss+density*width ss-density(end:-1:1)*width], ...
            [value value(end:-1:1)], colors{s},...
            'EdgeColor', colors{s});

        % plot upper and lower CI
        quartiles = [curr_CI(1),median(curr_error),curr_CI(2)];
        meanDensityWidth = interp1(value, density, quartiles(2))*width;
        hold on
        fill(ss+[-1,1,1,-1]*.02, ...
            [quartiles(1) quartiles(1) quartiles(3) quartiles(3)],...
            [0 0 0], 'EdgeColor', [0 0 0]);
        scatter(ss,quartiles(2),30,[1 1 1],'filled','black') % plot mean
        scatter(ss,quartiles(2),20,[1 1 1],'filled','white') % plot mean
    end

    plot(linspace(0.5,2.5,100),18*ones(1,100),'k:') % plot line at chance
    ylabel('Decoder error (cm)')
    ylim([0 36])
    xticklabels([{[]},stageList])
    title(categories{s})
    set(gca,'FontSize',14);


    set(gcf,'Position',[450 450 1000 400])
    set(gcf,'color',[1 1 1])
end