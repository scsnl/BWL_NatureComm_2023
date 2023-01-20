clc
clear all 
close all

Task = 'task';
Group = 'Chronos';
NumbROIs = 9;

%add relevant paths
addpath(fullfile(pwd,'/utils'));
% Load task condition 
load(fullfile(pwd, '/data/', '/stimulation_condition.mat'));
% Load data
load(fullfile(pwd, '/data/', sprintf('timeseries_ratAI_Opto_%s.mat', Group)));
% Load model
load(fullfile(pwd, '/model/', sprintf('BSDS_model_ratAI_Opto_%s.mat', Group)));
% Load covariance matrix
load(fullfile(pwd,'/model', sprintf('BSDS_covMtx_ratAI_Opto_%s.mat', Group)));


%% Data parameters
TR = 1; 
num.ROI = NumbROIs;
num.Subj = length(model.temporal_evolution_of_states); % number of subjects.
num.Vol = length(model.temporal_evolution_of_states{1,1}); % length of timeseries
num.State = length(unique(cell2mat(model.temporal_evolution_of_states)));
num.Run = 1;
num.Trials = 8; %number of ON stimulation blocks
ROI_names = {'PrL' 'AI', 'CG', 'RSC-2.9mm', 'RSC-3.9mm', 'RSC-4.9mm', 'RSC-5.9mm', 'RSC-6.9mm', 'RSC-7.8mm'};
ROI_names_excld_AI = {'PrL', 'CG', 'RSC-2.9mm', 'RSC-3.9mm', 'RSC-4.9mm', 'RSC-5.9mm', 'RSC-6.9mm', 'RSC-7.8mm'};
Colorstrings = {'#D62828','#E5E5E5','#003049', '#CFCFCF','#FFAC0B'};
Colormap = [0 48 73; 214 40 40; 255 172 11; 207 207 207; 229 229 229]./255;
reorder=[3,1,5,4,2]; %reorder the states for the visualization purpose.
ONblock.begin=[11, 111, 211, 311, 411, 511, 611, 711]; %indexes where ON blocks begin and end.
ONblock.end=[30, 130, 230, 330, 430, 530, 630, 730];

%% Check timeseries data of AI
for subj=1:num.Subj
    AI_ts(subj,:)=data{1,subj}(2,:);
end
figure;
subplot(2,1,1); plot(stimulation_cond);
xlim([0 length(stimulation_cond)]);
title('Stimulation protocol');
subplot(2,1,2); 
plot(mean(AI_ts));
xlim([0 length(stimulation_cond)]);
title('Time course of AI averaged across 9 Chronos rats');

%% Reshape group stats (i.e., temporal evolution of state, posterior probability)
% reshape temporal evolution of state into [subject x run] x [time] & identify dominant states
grpTempEvol = reshape(cell2mat(model.temporal_evolution_of_states),[num.Vol,num.Subj*num.Run])';
dominantStates = unique(grpTempEvol(:));
num.State=length(unique(cell2mat(model.temporal_evolution_of_states)));

% relabel states for ease of computation
grpTempEvol_relabel = zeros(size(grpTempEvol));
for relabel = 1:num.State
    grpTempEvol_relabel(grpTempEvol == dominantStates(relabel)) = relabel;
end

% reshape subjects' posterior probability(pp)
subjStatePP = zeros(num.Subj*num.Run, num.Vol, num.State); % each subject's posterior prob.
idx = 1;
for subj=1:num.Subj
    for run = 1:num.Run
        subjStatePP(idx,:,:) = model.posterior_probabilities{1,subj}((num.Vol*run-(num.Vol-1):num.Vol*run), dominantStates);
        idx = idx+1;
    end
end
grpStatePP = squeeze(mean(subjStatePP,1));% compute group posterior probability by averaging subjects posterior prob.


%% Figure3a: Plot task design & time-varying posterior probability & temporal evolution of the Chronos rats
figure;
subplot(3,1,1);
plot(stimulation_cond); xlim([0, num.Vol]); 
title('Task design');

subplot(3,1,2);
plot(grpStatePP);
for i = 1:num.State
    hold on; plot(grpStatePP(:,i), 'Color', Colorstrings{i}, 'LineWidth', 2);
end
ylim([0,1]); xlim([0, num.Vol]); xlabel('Time(s)', 'FontSize', 12); ylabel('posterior probability', 'FontSize', 12);
title('Posterior probability of latent brain states');

subplot(3,1,3);
hAxes = gca; imagesc(hAxes, grpTempEvol_relabel); colormap( hAxes, Colormap); 
title('Temporal Evolution of latent brain states');

%% Occupancy rates in OFF & ON periods
ON_periods = find(stimulation_cond==1);
OFF_periods = find(stimulation_cond==0);

OR_ON = zeros(num.Subj,num.State); OR_OFF = zeros(num.Subj,num.State); 

for subj=1:num.Subj
    tempEvol_ON=grpTempEvol_relabel(subj,ON_periods);
    [OR_ON(subj,:), ~, ~]=summary_stats_fast(tempEvol_ON, 1:num.State);
    tempEvol_OFF=grpTempEvol_relabel(subj,OFF_periods);
    [OR_OFF(subj,:), ~, ~]=summary_stats_fast(tempEvol_OFF, 1:num.State);
end
% Determine ON-State, OFF-State, and Transition-State
[~,ON_State]=max(mean(OR_ON)); [~,OFF_State]=max(mean(OR_OFF));
OR=[1:num.State; mean(OR_OFF)]'; OR([ON_State, OFF_State],:)=[]; OR_order=sortrows(OR, 2, 'descend');
Trans_State=OR_order(1,1);


%% Figure 3d: Plot occupancy rates of the Chronos rats
OR_ON = OR_ON(:, reorder); %reorder the states for the visualization purpose.
OR_OFF = OR_OFF(:, reorder); %reorder the states for the visualization purpose.
 
figure;
subplot(1,2,1);
b=bar(1:num.State, mean(OR_OFF)*100, 'FaceColor', 'flat');
for k=1:num.State; b.CData(k,:) = Colormap(k,:); end
hold on; errorbar(1:num.State, mean(OR_OFF).*100, (std(OR_OFF)/sqrt(num.Subj)).*100, 'LineStyle', 'none', 'Color', 'k');  % plot standard error.
set(gca, 'xticklabel', {'S1', 'S2', 'S3', 'S4', 'S5'}); ylim([0 100]); ylabel('Percentage');
title('Occupancy rate during OFF condition');

subplot(1,2,2);
b=bar(1:num.State, mean(OR_ON)*100, 'FaceColor', 'flat');
for k=1:num.State; b.CData(k,:) = Colormap(k,:); end
hold on; errorbar(1:num.State, mean(OR_ON).*100, (std(OR_ON)/sqrt(num.Subj)).*100, 'LineStyle', 'none', 'Color', 'k'); % plot standard error.
set(gca, 'xticklabel', {'S1', 'S2', 'S3', 'S4', 'S5'}); ylim([0 100]); ylabel('Percentage');
title('Occupancy rate during ON condition');


%% Dynamic state transition at stimulation boundaries
%% Figure 5a: State switching probability at OFF to ON stimulation boundary
%--------Check stimulation boundaries: OFF to ON-------%
count=1;
subjStatSwitchProb = zeros(num.Subj*num.Trials, num.State, num.State); % each subject's state switching matrix
for subj=1:num.Subj
    for trial=1:num.Trials
        OFF2ON=[ONblock.begin(trial)-10 : ONblock.begin(trial)+10];
        % y-axis:state(t), x-axis:state(t+1)
        SwitchingMatrix = full(sparse(grpTempEvol_relabel(subj, OFF2ON(1):OFF2ON(end)-1),grpTempEvol_relabel(subj, OFF2ON(2):OFF2ON(end)),1,num.State,num.State));
        sum_SwitchingMatrix = sum(SwitchingMatrix,2);
        % Normalize the transition matrix
        for state=1:num.State
            SwitchingMatrix(state, :) = SwitchingMatrix(state,:)/sum_SwitchingMatrix(state);
        end
        subjStatSwitchProb(count,:,:) = SwitchingMatrix;
        [OR_OFF2ON(count,:),~,~]=summary_stats_fast(grpTempEvol_relabel(subj, OFF2ON), 1:num.State);
        count=count+1;
    end
end
OFF2ONSwitchProb = squeeze(nanmean(subjStatSwitchProb(1:num.Subj,:,:)));
OFF2ONSwitchProb = OFF2ONSwitchProb(reorder, reorder);

figure;
subplot(1,2,1);
imagesc(OFF2ONSwitchProb);
ylabel('state(t)'); xlabel('state(t+1)');
title('State transition probability');
caxis([0 1]);
textStrings = num2str(OFF2ONSwitchProb(:), '%0.2f');       % Create strings from the matrix values
textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
[x, y] = meshgrid(1:num.State);  % Create x and y coordinates for the strings
hStrings = text(x(:), y(:), textStrings(:), ...  % Plot the strings
                'HorizontalAlignment', 'center');
midValue = mean(get(gca, 'CLim'));  % Get the middle value of the color range
% Choose white or black for the text color of the strings so they can be easily seen over the background color
textColors = repmat(OFF2ONSwitchProb(:) > midValue, 1, 3);  
set(hStrings, {'Color'}, num2cell(textColors, 2));  % Change the text colors
set(gca, 'XTick', 1:num.State, ...                             % Change the axes tick marks
         'XTickLabel', {'S1', 'S2', 'S3', 'S4', 'S5'}, ...  %   and tick labels
         'YTick', 1:num.State, ...
         'YTickLabel', {'S1', 'S2', 'S3', 'S4', 'S5'}, ...
         'TickLength', [0 0]);
mycolormap = customcolormap([0 0.5 1], {'#9d0142','#f66e45','#ffffbb'});
colorbar('southoutside');
colormap(mycolormap);

%% Figure 5a: State switching probability at ON to OFF stimulation boundary
%--------Check stimulation boundaries: ON to OFF-------%
count=1;
subjStatSwitchProb = zeros(num.Subj*num.Trials, num.State, num.State); % each subject's state switching matrix
for subj=1:num.Subj
    for trial=1:num.Trials
        ON2OFF=[ONblock.begin(trial)+10 : ONblock.end(trial)+10];
        % y-axis:state(t), x-axis:state(t+1)
        SwitchingMatrix = full(sparse(grpTempEvol_relabel(subj, ON2OFF(1):ON2OFF(end)-1),grpTempEvol_relabel(subj, ON2OFF(2):ON2OFF(end)),1,num.State,num.State));
        sum_SwitchingMatrix = sum(SwitchingMatrix,2);
        % Normalize the transition matrix
        for state=1:num.State
            SwitchingMatrix(state, :) = SwitchingMatrix(state,:)/sum_SwitchingMatrix(state);
        end
        subjStatSwitchProb(count,:,:) = SwitchingMatrix;
        [OR_ON2OFF(count,:),~,~]=summary_stats_fast(grpTempEvol_relabel(subj, ON2OFF), 1:num.State);
        count=count+1;
    end
end
ON2OFFSwitchProb = squeeze(nanmean(subjStatSwitchProb(1:num.Subj,:,:)));
ON2OFFSwitchProb = ON2OFFSwitchProb(reorder, reorder);

subplot(1,2,2);
imagesc(ON2OFFSwitchProb);
ylabel('state(t)'); xlabel('state(t+1)');
title('State transition probability during ON2OFF transition');
caxis([0 1]);
textStrings = num2str(ON2OFFSwitchProb(:), '%0.2f');       % Create strings from the matrix values
textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
[x, y] = meshgrid(1:num.State);  % Create x and y coordinates for the strings
hStrings = text(x(:), y(:), textStrings(:), ...  % Plot the strings
                'HorizontalAlignment', 'center');
midValue = mean(get(gca, 'CLim'));  % Get the middle value of the color range
% Choose white or black for the text color of the strings so they can be easily seen over the background color
textColors = repmat(ON2OFFSwitchProb(:) > midValue, 1, 3);  
set(hStrings, {'Color'}, num2cell(textColors, 2));  % Change the text colors
set(gca, 'XTick', 1:num.State, ...                             % Change the axes tick marks
         'XTickLabel', {'S1', 'S2', 'S3', 'S4', 'S5'}, ...  %   and tick labels
         'YTick', 1:num.State, ...
         'YTickLabel', {'S1', 'S2', 'S3', 'S4', 'S5'}, ...
         'TickLength', [0 0]);
mycolormap = customcolormap([0 0.5 1], {'#9d0142','#f66e45','#ffffbb'});
colorbar('southoutside');
colormap(mycolormap);


%% Activation level of each region in OFF, ON, and Transition States
ActivationLevel_ON=zeros(num.Subj, num.ROI);
ActivationLevel_OFF=zeros(num.Subj, num.ROI);
ActivationLevel_Trans=zeros(num.Subj, num.ROI);

for subj=1:num.Subj
    ActivationLevel_ON(subj,:)=model_subjwise.estimated_mean{1,subj}{1, dominantStates(ON_State)}';
    ActivationLevel_OFF(subj,:)=model_subjwise.estimated_mean{1,subj}{1, dominantStates(OFF_State)}';
    ActivationLevel_Trans(subj,:)=model_subjwise.estimated_mean{1,subj}{1, dominantStates(Trans_State)}';
end

[sig_t_mtx_ONvsOFF, p_mtx]=paired_ttest_activation_level((ActivationLevel_ON), (ActivationLevel_OFF), 0.01);
[sig_t_mtx_ONvsTrans, p_mtx]=paired_ttest_activation_level((ActivationLevel_ON), (ActivationLevel_Trans), 0.01);
[sig_t_mtx_TransvsOFF, p_mtx]=paired_ttest_activation_level((ActivationLevel_Trans), (ActivationLevel_OFF), 0.01);

%% Figure 4a and 6a,d
figure;
subplot(1,3,1);
bar([1:num.ROI], mean(ActivationLevel_OFF));
hold on;
errorbar([1:num.ROI], mean(ActivationLevel_OFF), (std(ActivationLevel_OFF)/sqrt(length(ActivationLevel_OFF))), 'LineStyle', 'none', 'Color', 'k'); % plot standard error.
set(gca,'xtick',[1:num.ROI],'xticklabel',ROI_names)
title('Activation level in OFF State');
ylim([0 1]);

subplot(1,3,2);
bar([1:num.ROI], mean(ActivationLevel_Trans));
hold on;
errorbar([1:num.ROI], mean(ActivationLevel_Trans), (std(ActivationLevel_Trans)/sqrt(length(ActivationLevel_Trans))), 'LineStyle', 'none', 'Color', 'k'); % plot standard error.
set(gca,'xtick',[1:num.ROI],'xticklabel',ROI_names)
title('Activation level in Trans State');
ylim([0 1]);

subplot(1,3,3);
bar([1:num.ROI], mean(ActivationLevel_ON));
hold on;
errorbar([1:num.ROI], mean(ActivationLevel_ON), (std(ActivationLevel_ON)/sqrt(length(ActivationLevel_ON))), 'LineStyle', 'none', 'Color', 'k'); % plot standard error.
set(gca,'xtick',[1:num.ROI],'xticklabel',ROI_names)
title('Activation level in ON State');
ylim([0 1]);

%% Figure 4b and 6b,e
figure;
subplot(1,3,1);
bar([1:num.ROI], sig_t_mtx_ONvsOFF);
set(gca,'xtick',[1:num.ROI],'xticklabel',ROI_names)
title('Activation level comparison: ON vs OFF States');
subplot(1,3,2);
bar([1:num.ROI], sig_t_mtx_ONvsTrans);
set(gca,'xtick',[1:num.ROI],'xticklabel',ROI_names)
title('Activation level comparison: ON vs Trans States');
subplot(1,3,3);
bar([1:num.ROI], sig_t_mtx_TransvsOFF);
set(gca,'xtick',[1:num.ROI],'xticklabel',ROI_names)
title('Activation level comparison: Trans vs OFF States');


%% Connectivity analysis: 

cov_subj=zeros(num.ROI, num.ROI, num.State, num.Subj);
partial_corr_subj=zeros(num.ROI, num.ROI, num.State, num.Subj);
partial_corr_subj_z=zeros(num.ROI, num.ROI, num.State, num.Subj);
for j= 1:length(dominantStates)
    grp_avg_fc=zeros(num.ROI, num.ROI);
    k = dominantStates(j);
    for subj=1:num.Subj
        est_cov = model_subjwise.estimated_covariance{1,subj}{1,k};
        cov_subj(:,:,j,subj) = est_cov;
        %Partial Correlation
        inv_est_cov = inv(est_cov);
        invD = inv(diag(sqrt(diag(inv_est_cov))));
        partial_corr_subj(:,:,j,subj) = -invD*inv_est_cov*invD;
        partial_corr_subj_z(:,:,j,subj)=log((1+partial_corr_subj(:,:,j,subj))./(1-partial_corr_subj(:,:,j,subj)));        
    end
end
Corr_OFF_z=squeeze(partial_corr_subj_z(:,:,OFF_State,:));
Corr_ON_z=squeeze(partial_corr_subj_z(:,:,ON_State,:));

pval=0.01;
[sig_t_mtx_ONvsOFF, p_mtx] = paired_ttest_corr_mtx(partial_corr_subj_z(:,:,ON_State,:), partial_corr_subj_z(:,:,OFF_State,:), pval);
[sig_t_mtx_ONvsTrans, p_mtx] = paired_ttest_corr_mtx(partial_corr_subj_z(:,:,ON_State,:), partial_corr_subj_z(:,:,Trans_State,:), pval);
[sig_t_mtx_TransvsOFF, p_mtx] = paired_ttest_corr_mtx(partial_corr_subj_z(:,:,Trans_State,:), partial_corr_subj_z(:,:,OFF_State,:), pval);

%% Figure 4c, 6c,f
figure;
subplot(1,3,1);
plot_connectivity_diff(sig_t_mtx_ONvsOFF, num.ROI, ROI_names,[-6 6], 'Conn of ON compared to OFF');
subplot(1,3,2);
plot_connectivity_diff(sig_t_mtx_ONvsTrans, num.ROI, ROI_names,[-6 6], 'Conn of ON compared to Trans');
subplot(1,3,3);
plot_connectivity_diff(sig_t_mtx_TransvsOFF, num.ROI, ROI_names,[-6 6], 'Conn of Trans compared to OFF');

%% Figure 4d: Comparing AI-other ROIs connection between ON & OFF states.
AI_RSCs_OFF=(squeeze(Corr_OFF_z(2,[1,3:NumbROIs],:))');
AI_RSCs_ON=(squeeze(Corr_ON_z(2,[1,3:NumbROIs],:))');
xaxis_range=[1:NumbROIs-1];
figure;
plot(xaxis_range, mean(AI_RSCs_OFF), '-ok');
hold on;
errorbar(xaxis_range, mean(AI_RSCs_OFF), (std(AI_RSCs_OFF)/sqrt(length(AI_RSCs_OFF))), 'LineStyle', 'none', 'Color', 'k'); % plot standard error.
set(gca,'xtick',xaxis_range,'xticklabel',ROI_names_excld_AI)
hold on;
plot(xaxis_range, mean(AI_RSCs_ON), '-or');
hold on;
errorbar(xaxis_range, mean(AI_RSCs_ON), (std(AI_RSCs_ON)/sqrt(length(AI_RSCs_ON))), 'LineStyle', 'none', 'Color', 'r'); % plot standard error.
set(gca,'xtick',xaxis_range,'xticklabel',ROI_names_excld_AI)
title('connection between AI and other ROIs');
ylim([-0.4, 0.4]);

%% Figure 4e: Comparing PrL-other ROIs connection between ON & OFF states.
ROI_names={'PrL' 'AI', 'CG', 'RSC-2.9mm', 'RSC-3.9mm', 'RSC-4.9mm', 'RSC-5.9mm', 'RSC-6.9mm', 'RSC-7.8mm'};
ROI_names_excld_PrL={'AI', 'CG', 'RSC-2.9mm', 'RSC-3.9mm', 'RSC-4.9mm', 'RSC-5.9mm', 'RSC-6.9mm', 'RSC-7.8mm'};

PrL_OFF=(squeeze(Corr_OFF_z(1,[2:NumbROIs],:))');
PrL_ON=(squeeze(Corr_ON_z(1,[2:NumbROIs],:))');
xaxis_range=[1:NumbROIs-1];
figure;
plot(xaxis_range, mean(PrL_OFF), '-ok');
hold on;
errorbar(xaxis_range, mean(PrL_OFF), (std(PrL_OFF)/sqrt(length(PrL_OFF))), 'LineStyle', 'none', 'Color', 'k'); % plot standard error.
set(gca,'xtick',xaxis_range,'xticklabel',ROI_names_excld_PrL)
hold on;
plot(xaxis_range, mean(PrL_ON), '-or');
hold on;
errorbar(xaxis_range, mean(PrL_ON), (std(PrL_ON)/sqrt(length(PrL_ON))), 'LineStyle', 'none', 'Color', 'r'); % plot standard error.
set(gca,'xtick',xaxis_range,'xticklabel',ROI_names_excld_PrL)
title('connection between PrL and other ROIs');
ylim([-0.4, 0.7]);

%% Figure 4f: Comparing RSC-6.8-other ROIs connection between ON & OFF states.
ROI_names={'PrL' 'AI', 'CG', 'RSC-2.9mm', 'RSC-3.9mm', 'RSC-4.9mm', 'RSC-5.9mm', 'RSC-6.9mm', 'RSC-7.8mm'};
ROI_names_excld_RSC={'PrL', 'AI', 'CG', 'RSC-2.9mm', 'RSC-3.9mm', 'RSC-4.9mm', 'RSC-5.9mm', 'RSC-7.8mm'};

RSC_OFF=(squeeze(Corr_OFF_z(8,[1:7,NumbROIs],:))');
RSC_ON=(squeeze(Corr_ON_z(8,[1:7,NumbROIs],:))');
xaxis_range=[1:NumbROIs-1];
figure;
plot(xaxis_range, mean(RSC_OFF), '-ok');
hold on;
errorbar(xaxis_range, mean(RSC_OFF), (std(RSC_OFF)/sqrt(length(RSC_OFF))), 'LineStyle', 'none', 'Color', 'k'); % plot standard error.
set(gca,'xtick',xaxis_range,'xticklabel',ROI_names_excld_RSC)
hold on;
plot(xaxis_range, mean(RSC_ON), '-or');
hold on;
errorbar(xaxis_range, mean(RSC_ON), (std(RSC_ON)/sqrt(length(RSC_ON))), 'LineStyle', 'none', 'Color', 'r'); % plot standard error.
set(gca,'xtick',xaxis_range,'xticklabel',ROI_names_excld_RSC)
title('connection between RSC-6.8 and other ROIs');
ylim([-0.2, 1]);








