function [plotResults] = PC_Position_Decoder_allMice(plotFlag, fname)
%%% Calculate decoder error across stages and mice
%   Randomly sample using minimum number of trials (across all stages/mice)
%   Randomly sample using minimum number of cells (across stages)
%   Written: MG 231005, Latest update: NW 240611

clearvars -except plotFlag fname
close all

if nargin==0 || isempty(plotFlag)
    plotFlag=1;
end
if nargin < 2 || isempty(fname)
    %% load data
    disp('Select place cell responses...')
    [lname, lpath] = uigetfile('*.mat');
    fname = fullfile(lpath,lname);
end

%% load data
% file should be structured as an S x M cell array, where S is estrous stages and M is mouse.
% each cell should be a 1x3 cell array, where 3 is environments A > B > A'
% Within each environment is an LxBxC array, where L is # laps, B is # bins, and C is # cells
laps_by_animal = importdata(fname);

%% user parameters
smoothFactor = 1;       % smooth maximum likelihood vector to improve position estimates (default: 1 = no smoothing)
numMouse = 6;           % number of mice used in experiment (default: n = 6)
numStage = 4;           % number of mice used in experiment (default: 4 stages)
numTrial = 7;           % min number of trials used in experiment default: 10 trials)

%% extract matrix with number of cells per mouse/stage, calculate minimum number of cells to use
numCelMat = zeros(numStage,numMouse);
for s = 1:numStage
    for m = 1:numMouse
        numCelMat(s,m) = size(laps_by_animal{s,m}{1},3);
    end
end
numPos = size(laps_by_animal{s,m}{1},2);
numCel = min(numCelMat);        % mininum number of cells across stages per mouse (to keep stages equal)

%% initialize likelihood and error matrices

likelihoodMat_B = zeros(numPos,numPos,numStage);
likelihoodMat_Ap = zeros(numPos,numPos,numStage);

error_B = zeros(numTrial*numPos,numStage);
error_B_avg = zeros(numPos,numStage);

error_Ap = zeros(numTrial*numPos,numStage);
error_Ap_avg = zeros(numPos,numStage);

%% Decoding analysis
for s = 1:numStage    % run decoding across all four stages
    
    resp_mat_A = [];
    resp_mat_B = [];
    resp_mat_Ap = [];
    
    %% populate response matrices
    for m = 1:numMouse
        
        curr_session = laps_by_animal{s,m};
        
        % randomize selected cells
        rand_cellvec = randperm(size(curr_session{1},3));
        rand_cell_idx = rand_cellvec(1:numCel(m));
        
        % randomize trials and add to response matrix
        rand_trialvec = randperm(size(curr_session{1},1));
        rand_trial_idx = rand_trialvec(1:numTrial);
        resp_mat_A = cat(3,resp_mat_A, curr_session{1}(rand_trial_idx,:,rand_cell_idx));
        
        rand_trialvec = randperm(size(curr_session{2},1));
        rand_trial_idx = rand_trialvec(1:numTrial);
        resp_mat_B = cat(3,resp_mat_B, curr_session{2}(rand_trial_idx,:,rand_cell_idx));
        
        rand_trialvec = randperm(size(curr_session{3},1));
        rand_trial_idx = rand_trialvec(1:numTrial);
        resp_mat_Ap = cat(3,resp_mat_Ap, curr_session{3}(rand_trial_idx,:,rand_cell_idx));
        
    end
    
    % normalize responses
    for c = 1:sum(numCel)
        
        minval = min([min(min(resp_mat_A(:,:,c))) min(min(resp_mat_B(:,:,c))) min(min(resp_mat_Ap(:,:,c)))]);
        maxval = max([max(max(resp_mat_A(:,:,c))) max(max(resp_mat_B(:,:,c))) max(max(resp_mat_Ap(:,:,c)))]);
        
        resp_mat_A(:,:,c) = resp_mat_A(:,:,c)-minval;
        resp_mat_A(:,:,c) = resp_mat_A(:,:,c)/maxval;
        
        resp_mat_B(:,:,c) = resp_mat_B(:,:,c)-minval;
        resp_mat_B(:,:,c) = resp_mat_B(:,:,c)/maxval;
        
        resp_mat_Ap(:,:,c) = resp_mat_Ap(:,:,c)-minval;
        resp_mat_Ap(:,:,c) = resp_mat_Ap(:,:,c)/maxval;
        
    end
    
    
    
    %% initialize matrices
    actualPos_B = zeros(1,numTrial*numPos);
    estimatedPos_B = zeros(1,numTrial*numPos);
    
    actualPos_Ap = zeros(1,numTrial*numPos);
    estimatedPos_Ap = zeros(1,numTrial*numPos);
    
    %% Extract models and responses for decoding (remove NaNs)
    envA_resp = resp_mat_A;
    envA_resp(isnan(envA_resp)) = 0;
    modelResp = squeeze(mean(envA_resp,1))';
    
    envB_resp = resp_mat_B;
    envB_resp(isnan(envB_resp)) = 0;
    decoderResp_B = envB_resp;
    
    envAp_resp = resp_mat_Ap;
    envAp_resp(isnan(envAp_resp)) = 0;
    decoderResp_Ap = envAp_resp;
    
    %% A>B single trial decoding
    for lap = 1:numTrial
        for pos = 1:numPos
            
            popVec_B = squeeze(decoderResp_B(lap,pos,:));
            
            likelihoodVec = zeros(1,72);
            
            for c = 1:sum(numCel)
                likelihoodVec = likelihoodVec+(popVec_B(c)*modelResp(c,:));
            end
            
            likelihoodVec = smooth(likelihoodVec,smoothFactor);
            [~,MaxProbPos] = max(likelihoodVec);
            
            actualPos_B((lap-1)*numPos+pos) = pos;
            estimatedPos_B((lap-1)*numPos+pos) = MaxProbPos;
            err_diff1 = max(MaxProbPos,pos)-min(MaxProbPos,pos);
            err_diff2 = min(MaxProbPos,pos)+72-max(MaxProbPos,pos);
            error_B((lap-1)*numPos+pos,s) = min(err_diff1,err_diff2);
        end
    end
    
    %% A>A' single trial decoding
    for lap = 1:numTrial
        for pos = 1:numPos
            
            popVec_Ap = squeeze(decoderResp_Ap(lap,pos,:));
            
            likelihoodVec = zeros(1,72);
            
            for c = 1:sum(numCel)
                likelihoodVec = likelihoodVec+(popVec_Ap(c)*modelResp(c,:));
            end
            
            likelihoodVec = smooth(likelihoodVec,smoothFactor);
            [~,MaxProbPos] = max(likelihoodVec);
            
            actualPos_Ap((lap-1)*numPos+pos) = pos;
            estimatedPos_Ap((lap-1)*numPos+pos) = MaxProbPos;
            err_diff1 = max(MaxProbPos,pos)-min(MaxProbPos,pos);
            err_diff2 = min(MaxProbPos,pos)+72-max(MaxProbPos,pos);
            error_Ap((lap-1)*numPos+pos,s) = min(err_diff1,err_diff2);
        end
    end
    
    %% A>B average trial decoding
    for pos = 1:numPos
        
        popVec_B = squeeze(median(decoderResp_B(:,pos,:),1));
        
        likelihoodVec = zeros(1,72);
        
        for c = 1:sum(numCel)
            likelihoodVec = likelihoodVec+(popVec_B(c)*modelResp(c,:));
        end
        
        likelihoodVec = smooth(likelihoodVec,smoothFactor);
        likelihoodMat_B(:,pos,s) = likelihoodVec'/mean(likelihoodVec);
        [~,MaxProbPos] = max(likelihoodVec);
        
        err_diff1 = max(MaxProbPos,pos)-min(MaxProbPos,pos);
        err_diff2 = min(MaxProbPos,pos)+72-max(MaxProbPos,pos);
        error_B_avg(pos,s) = min(err_diff1,err_diff2);
    end
    
    %% A>A' average trial decoding
    for pos = 1:numPos
        
        popVec_Ap = squeeze(median(decoderResp_Ap(:,pos,:),1));
        
        likelihoodVec = zeros(1,72);
        
        for c = 1:sum(numCel)
            likelihoodVec = likelihoodVec+(popVec_Ap(c)*modelResp(c,:));
        end
        
        likelihoodVec = smooth(likelihoodVec,smoothFactor);
        likelihoodMat_Ap(:,pos,s) = likelihoodVec'/mean(likelihoodVec);
        [~,MaxProbPos] = max(likelihoodVec);
        
        err_diff1 = max(MaxProbPos,pos)-min(MaxProbPos,pos);
        err_diff2 = min(MaxProbPos,pos)+72-max(MaxProbPos,pos);
        error_Ap_avg(pos,s) = min(err_diff1,err_diff2);
    end
    
    %% Plot population decoding
    
    if plotFlag==1
        
        %% Plot "single" trial decoding
        
        figure(2*(s-1)+1)
        
        subplot(7,1,1:2)
        plot(actualPos_B,'k:')
        hold on
        scatter(1:length(estimatedPos_B),estimatedPos_B,'bo','filled')
        xlim([1 length(estimatedPos_B)])
        ylim([0 numPos+1])
        title('Prediction of Environment B position using Environment A model')
        
        subplot(7,1,3)
        plot(error_B(:,s),'r')
        xlim([1 length(estimatedPos_B)])
        ylim([0 numPos/2+1])
        title('Prediction error (Environment B)')
        
        subplot(7,1,5:6)
        plot(actualPos_Ap,'k:')
        hold on
        scatter(1:length(estimatedPos_Ap),estimatedPos_Ap,'bo','filled')
        xlim([1 length(estimatedPos_Ap)])
        ylim([0 numPos+1])
        title('Prediction of Environment A prime position using Environment A model')
        
        subplot(7,1,7)
        plot(error_Ap(:,s),'r')
        xlim([1 length(estimatedPos_Ap)])
        ylim([0 numPos/2+1])
        title('Prediction error (Environment A prime)')
        
        set(gcf,'color',[1 1 1])
        set(gcf,'position',[500 100 1000 600])
        
        %% Plot average trial decoding
        
        figure(2*s)
        
        subplot(3,3,[1 4])
        imagesc(likelihoodMat_B(:,:,s))
        set(gca,'YDir','normal')
        hold on
        plot(1:numPos,'k:')
        xlabel('Actual position')
        ylabel('Estimate position (maximum likelihood)')
        axis square
        title('Prediction of Environment B position using Environment A model')
        
        subplot(3,3,7)
        plot(error_B_avg(:,s),'r')
        hold on
        plot(1:numPos,18*ones(1,numPos),'k:')
        xlim([1 numPos])
        ylim([0 numPos/2+1])
        xlabel('Position (cm)')
        ylabel('Decoder error (cm)')
        title('Prediction error (Environment B)')
        
        subplot(3,3,[2 5])
        imagesc(likelihoodMat_Ap(:,:,s))
        set(gca,'YDir','normal')
        hold on
        plot(1:numPos,'k:')
        xlabel('Actual position')
        ylabel('Estimate position (maximum likelihood)')
        axis square
        title('Prediction of Environment A prime position using Environment A model')
        
        subplot(3,3,8)
        plot(error_Ap_avg(:,s),'r')
        hold on
        plot(1:numPos,18*ones(1,numPos),'k:')
        xlim([1 numPos])
        ylim([0 numPos/2+1])
        xlabel('Position (cm)')
        ylabel('Decoder error (cm)')
        title('Prediction error (Environment A prime)')
        
        subplot(3,3,[3 6])
        scatter(ones(1,numPos)+(rand(1,numPos)-0.5)/4,error_B_avg(:,s),25,[0 0.75 0.75])
        hold on
        scatter(1,mean(error_B_avg(:,s)),50,[0 0.5 0.5],'filled')
        negSE = mean(error_B_avg(:,s))-std(error_B_avg(:,s))/sqrt(numPos);
        posSE = mean(error_B_avg(:,s))+std(error_B_avg(:,s))/sqrt(numPos);
        plot(ones(1,100),linspace(negSE,posSE,100),'linewidth',2,'color',[0 0.5 0.5])
        scatter(2*ones(1,numPos)+(rand(1,numPos)-0.5)/4,error_Ap_avg(:,s),25,[0.75 0 0.75])
        scatter(2,mean(error_Ap_avg(:,s)),50,[0.5 0 0.5],'filled')
        negSE = mean(error_Ap_avg(:,s))-std(error_Ap_avg(:,s))/sqrt(numPos);
        posSE = mean(error_Ap_avg(:,s))+std(error_Ap_avg(:,s))/sqrt(numPos);
        plot(2*ones(1,100),linspace(negSE,posSE,100),'linewidth',2,'color',[0.5 0 0.5])
        plot(linspace(0.5,2.5,100),18*ones(1,100),'k:')
        ylabel('Decoder error (cm)')
        xticks([0.5 1 1.5 2 2.5])
        xticklabels({' ','A>B',' ','A>Ap',' '})
        title('Decoder error')
        
        set(gcf,'color',[1 1 1])
        set(gcf,'position',[600 300 1000 450])
        
        %% print mean errors
        
        stageList = {'Diestrus','Proestrus','Estrus','Metestrus'};
        disp([stageList{s} ':'])
        disp(['A>B error (average trial): ' num2str(mean(error_B_avg(:,s))) '+/-' num2str(std(error_B_avg(:,s))) ' mean +/- s.d.'])
        disp(['A>Ap error (average trial): ' num2str(mean(error_Ap_avg(:,s))) '+/-' num2str(std(error_Ap_avg(:,s))) ' mean +/- s.d.'])
        disp(' ')
        
    end
    
end

%% Save to structure for iterations
plotResults.numPos = numPos;
plotResults.error_B = error_B_avg;
plotResults.error_Ap = error_Ap_avg;
plotResults.likelihood_B = likelihoodMat_B;
plotResults.likelihood_Ap = likelihoodMat_Ap;